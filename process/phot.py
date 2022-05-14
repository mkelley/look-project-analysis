#!/usr/bin/env python3
import os
from glob import glob
import argparse
import logging
from collections import defaultdict
from itertools import groupby
import warnings
from mskpy.image.core import rebin

import numpy as np
import scipy.ndimage as nd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import log as astropy_log
from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, FITSFixedWarning
import sep
from sbpy.data import Ephem, Orbit, TimeScaleWarning
from sbpy.activity import phase_HalleyMarcus

from calviacat import RefCat2
from mskpy import gcentroid
from mskpy.photometry.outbursts import CometaryTrends
from lco_util import *

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('files', nargs='*', help='measure these files')
parser.add_argument('-f', '--force', action='store_true',
                    help='re-measure and overwrite existing data')
parser.add_argument('--debug', action='store_true',
                    help='enable debugging messages for the screen')
parser.add_argument('-q', action='store_true', help='quiet mode, reduces messages on the screen')
args = parser.parse_args()

if not args.debug:
    warnings.simplefilter('ignore', FITSFixedWarning)
    warnings.simplefilter('ignore', RuntimeWarning)

logger = logging.getLogger('lco-phot')
if len(logger.handlers) == 0:
    handler = logging.StreamHandler()
    if args.q:
        handler.setLevel(logging.WARNING)
    elif args.debug:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    logger.addHandler(handler)

    handler = logging.FileHandler('phot.log')
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)

# prevents logging of fitter warnings from calviacat
astropy_log.setLevel('ERROR')

if not os.path.exists('zeropoints'):
    os.system('mkdir zeropoints')

if not os.path.exists('backgrounds'):
    os.system('mkdir backgrounds')

rho_arcsec = [2, 5, 10, 12, 20]
rho_km = [5e3, 1e4, 2e4]
rho_labels = (
    [str(r) for r in rho_arcsec] +
    ['{}k'.format(str(int(r // 1000))) for r in rho_km]
)
color_rho = '5'

# skip_objects = ['C/2021 C4']  # not in horizons

files = glob('e91/*/*fz')
if os.path.exists('phot-skip.list'):
    skip = ascii.read('phot-skip.list')
else:
    skip = Table(names=('file', 'reason'), dtype=('U64', 'U128'))


bgopts = dict(bw=256, bh=256, fw=5, fh=5)


def set_table_formats(tab):
    tab['fracday'].format = '{:.5f}'
    tab['tmtp'].format = '{:.3f}'
    tab['rh'].format = '{:.3f}'
    tab['delta'].format = '{:.3f}'
    tab['phase'].format = '{:.3f}'
    tab['zp'].format = '{:.4f}'
    tab['color cor'].format = '{:.4f}'
    tab['zp err'].format = '{:.4f}'
    tab['exptime'].format = '{:.1f}'
    tab['airmass'].format = '{:.3f}'
    tab['seeing'].format = '{:.2f}'
    tab['rho10k'].format = '{:.2f}'
    tab['cx'].format = '{:.1f}'
    tab['cy'].format = '{:.1f}'
    tab['dc'].format = '{:.1f}'
    tab['bgarea'].format = '{:.0f}'
    tab['bg'].format = '{:.5g}'
    tab['bgsig'].format = '{:.5g}'
    for r in rho_labels:
        tab['flux{}'.format(r)].format = '{:.5g}'
        tab['m{}'.format(r)].format = '{:.3f}'
        tab['merr{}'.format(r)].format = '{:.3f}'
    return tab


if os.path.exists('phot.txt') and not args.force:
    tab = set_table_formats(ascii.read('phot.txt'))
else:
    tab = {'file': []}

# ephemeris offset
offset = defaultdict(lambda: {
    'ra': 0 * u.arcsec,
    'dec': 0 * u.arcsec
})
offset['141P'] = {
    'ra':  1.4 * u.arcsec,
    'dec': -0.9 * u.arcsec
}

# centroiding box (box sizes are binned pixels)
centroid_options = defaultdict(lambda: {
    'binning': 1,
    'boxes': [17, 11]
})
centroid_options['C/2014 UN271'] = {
    'binning': 2,
    'boxes': None
}
# centroid_options['7P'] = {
#     'binning': 1,
#     'boxes': [9, 5]
# }
# centroid_boxes = defaultdict(lambda: [17, 11])
# centroid_boxes['2014 UN271'] = [21]

new_objects = set()

for f in files:
    basename = os.path.basename(f)[:-8]
    if f in tab['file'] and basename in skip['file']:
        # this file was added to the skip list after photometry; remove it
        if isinstance(tab, dict):
            i = tab['file'].index(f)
            del tab['file'][i]
        else:
            i = tab['file'] != f
            tab = tab[i]

    if basename in skip['file'] or f in tab['file']:
        continue

    bgf = 'backgrounds/{}.fits'.format(basename)

    try:
        hdu = fits.open(f)
    except OSError:
        logger.exception('Error opening %s.', basename)
        skip.add_row((basename, 'Error opening file'))
        skip.write('phot-skip.list',
                   format='ascii.csv', overwrite=True)
        continue
    h = hdu['sci'].header

    target = rename_target.get(h['OBJECT'], h['OBJECT'])
    # if target != '2014 UN271':
    #     continue

    bpm = hdu['bpm'].data != 0
    im = hdu['sci'].data + 0
    # do not use gain, already in e-
    err = np.sqrt((im + h['BIASLVL']) + h['RDNOISE']**2)

    if 'cat' not in hdu:
        logger.exception('%s missing photometry catalog.', basename)
        skip.add_row((basename, 'Missing photometry catalog'))
        skip.write('phot-skip.list', format='ascii.csv',
                   overwrite=True)
        continue
    
    phot = Table(hdu['cat'].data)
    if 'ra' not in phot.colnames:
        logger.error('Adding %s to skip file list: WCS probably missing.',
                     basename)
        skip.add_row((basename, 'WCS probably missing'))
        skip.write('phot-skip.list',
                   format='ascii.csv', overwrite=True)
        continue

    new_objects.add(target)
    ps = h['PIXSCALE'] * u.arcsec

    # first determine if target is in the image
    midtime = Time(h['DATE-OBS']) + h['EXPTIME'] / 2 * u.s
    opts = dict(epochs=midtime, cache=True, location=locations[basename[:3]])
    eph = None
    try:
        if target.startswith('C/') or target.startswith('P/') or target.endswith('P'):
            opts2 = dict(
                id_type='designation',
                closest_apparition=True,
                no_fragments=True,
            )
        else:
            opts2 = dict(
                id_type="smallbody"
            )

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', (UserWarning, TimeScaleWarning))
            eph = Ephem.from_horizons(target, **opts, **opts2)
            orb = Orbit.from_horizons(target, epochs=opts['epochs'],
                                      cache=opts['cache'], **opts2)

        eph['rh'] = [np.sign(row['rdot']) * row['rh'] for row in eph]
        Tp = orb['Tp'].jd
    except:
        pass

    try:
        if eph is None:
            eph = Ephem.from_mpc(target, **opts)
            Tp = np.ma.masked
    except:
        continue

    # target location
    wcs = WCS(h)
    c = wcs.all_world2pix(eph['RA'] + offset[target]['ra'],
                          eph['DEC'] + offset[target]['dec'], 0)
    if any(np.isnan(c).ravel()):
        logger.error('Adding %s to skip file list: WCS transformation error.'
                     '  Ephemeris: %s, %s', basename,
                     eph['RA'] + offset[target]['ra'],
                     eph['DEC'] + offset['dec'])
        skip.add_row((basename, 'WCS transformation error'))
        skip.write('phot-skip.list',
                   format='ascii.csv', overwrite=True)
        continue

    gx, gy = np.array(c).ravel()

    if min(gx, gy) < 0 or gx > im.shape[1] or gy > im.shape[0]:
        logger.error('Adding %s to skip file list: Target not in frame.',
                     basename)
        skip.add_row((basename, 'Target not in frame'))
        skip.write('phot-skip.list',
                   format='ascii.csv', overwrite=True)
        continue

    # source mask
    if not os.path.exists(bgf):
        # quick bg filter
        mask = bpm.copy()
        bkg = sep.Background(im, mask=mask, **bgopts)
        _im = im - bkg.back()

        mask += nd.grey_closing(_im, 5) > bkg.globalrms * 2
    else:
        mask = fits.getdata(bgf, 'mask').astype(bool)

    source_mask = mask
    del mask

    # calibrate
    lco = SkyCoord(phot['ra'], phot['dec'], unit='deg')

    rc2 = RefCat2('cat.db', min_matches=10, logger=logger)

    if len(rc2.search(lco)[0]) < 100:
        rc2.fetch_field(lco)

    try:
        objids, distances = rc2.xmatch(lco)
    except TypeError:
        logger.error('Adding %s to skip file list: Error in xmatch.',
                     basename)
        skip.add_row((basename, 'Error in catalog xmatch'))
        skip.write('phot-skip.list',
                   format='ascii.csv', overwrite=True)
        continue

    m_inst = -2.5 * np.log10(phot['flux'])
    m_err = phot['fluxerr'] / phot['flux'] * 1.0857
    catfilt = filters[h['FILTER']]

    # if h['FILTER'] == 'gp':
    #    C = -0.120
    # elif h['FILTER'] == 'rp':
    #    C = 0.029
    # else:
    #    C = None

    C = color_corrections.get(h['filter'])
    try:
        zp, C, zp_unc, m, gmr, gmi = rc2.cal_color(
            objids, m_inst, catfilt, 'g-r', C=C, mlim=[13, 18],
            gmi_lim=[0.2, 3.0])
    except Exception as e:
        logger.exception('Calibration error with %s', basename)
        skip.add_row((basename, 'Calibration error.'))
        skip.write('phot-skip.list',
                   format='ascii.csv', overwrite=True)
        continue
    n_cal = (~m.mask).sum()

    # plot
    fig = plt.figure(1)
    fig.clear()
    ax = fig.gca()

    ax.scatter(gmr, m - m_inst, marker='.', color='k')

    x = np.linspace(0, 1.5)
    label = '{:.4f} + {:.4f} ($g-r$)'.format(zp, C)
    ax.plot(x, C * x + zp, 'r-', label=label)

    plt.setp(ax, xlabel='$g-r$ (mag)', ylabel=r'$m-m_{\rm inst}$ (mag)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('zeropoints/{}.png'.format(basename))

    # 2D background
    if not os.path.exists(bgf):
        bkg = sep.Background(im, mask=source_mask, **bgopts)
        bghdu = fits.HDUList()
        bghdu.append(fits.PrimaryHDU(bkg.back().astype(np.float32)))
        bghdu.append(fits.ImageHDU(source_mask.astype(np.uint8),
                                   name='mask'))
        bghdu[0].header['median'] = (bkg.globalback, 'global median')
        bghdu[0].header['medrms'] = (bkg.globalrms, 'global rms')
        bghdu.writeto(bgf, overwrite=True)

    bgim, bgh = fits.getdata(bgf, header=True)
    im -= bgim
    bg = 0
    bgsig = bgh['medrms']
    bgarea = np.sum(~source_mask)

    binning = centroid_options[target]['binning']
    boxes = centroid_options[target]['boxes']
    if boxes is None:
        cy = gy
        cx = gx
        dc = np.ma.masked
    else:
        binned = rebin(im, -binning, trim=True)
        gyb = (gy + 0.5) / binning - 0.5
        gxb = (gx + 0.5) / binning - 0.5
        cyb, cxb = gyb, gxb
        for box in boxes:
            try:
                cyb, cxb = gcentroid(binned, (cyb, cxb), box=box)
                dc = np.hypot(cxb - gxb, cyb - gyb) * binning
            except:
                dc = np.ma.masked

        cy = (cyb + 0.5) * binning - 0.5
        cx = (cxb + 0.5) * binning - 0.5

    # for box in centroid_boxes.get(target):
    #     try:
    #         cy, cx = gcentroid(im, (cy, cx), box=box)
    #         dc = np.hypot(cx - gx, cy - gy)
    #     except:
    #         dc = np.ma.masked

    # background annulus
    # bgap = pu.CircularAnnulus((cx, cy), 30, 50)
    # mask = bgap.to_mask(method='center')
    # b = mask.multiply(im)
    # mms = sigma_clipped_stats(b[b != 0])
    # bg = mms[1]
    # #bg = 0
    # bgsig = mms[2]
    # bgarea = bgap.area

    rho10k = (1e4 * u.km / (725 * u.km / u.arcsec / u.au) / eph['delta']).to(
        u.arcsec).value[0]
    rho = [r * u.arcsec for r in rho_arcsec]
    rho.extend([r / 725 / eph['delta'].value * u.arcsec for r in rho_km])
    rho = u.Quantity(rho, u.arcsec)
    rap = (rho / ps).value
    area = np.pi * rap**2

    flux, err, flags = sep.sum_circle(im, [cx], [cy], rap, err=err)
    flux -= area * bg
    err2 = err**2 + area * bgsig**2 * (1 + area / bgarea)

    mc = -2.5 * np.log10(flux) + zp + C * assumed_gmr[target]
    mc_err = np.sqrt(err2 / (flux * 1.0857)**2 + zp_unc**2)

    t = Time(midtime)
    fracday = t.mjd % 1
    tmtp = t.jd - Tp

    row = ([f, target, midtime.iso[:-4], fracday, tmtp, eph['ra'].value[0],
            eph['dec'].value[0], eph['r'].value[0],
            eph['delta'].value[0], eph['alpha'].value[0], h['FILTER'], catfilt, n_cal,
            zp, C, zp_unc, h['EXPTIME'], h['AIRMASS'], h['L1FWHM'], rho10k,
            cx, cy, dc, bgarea, bg, bgsig] + list(flux) + list(mc) + list(mc_err))

    if isinstance(tab, dict):
        names = ['file', 'target', 'date', 'fracday', 'tmtp', 'ra', 'dec',
                 'rh', 'delta', 'phase', 'filter', 'catalog filter', 'N cal',
                 'zp', 'color cor', 'zp err', 'exptime', 'airmass', 'seeing',
                 'rho10k', 'cx', 'cy', 'dc', 'bgarea', 'bg', 'bgsig']
        names.extend(['flux{}'.format(r) for r in rho_labels])
        names.extend(['m{}'.format(r) for r in rho_labels])
        names.extend(['merr{}'.format(r) for r in rho_labels])
        tab = Table(rows=[row], names=names, masked=True)
        tab.meta['comments'] = [
            Time.now().iso,
            'Filters calibrated to PS1: ' + str(filters),
            'Color corrected assuming g-r = ',
        ] + ['  {}: {}'.format(k, v) for k, v in assumed_gmr.items()] + [
            'm2, m5, etc. use 2", 5", etc radius apertures',
            'm10k, etc. use 10,000 km etc. radius apertures'
        ]
        tab = set_table_formats(tab)
    else:
        tab.add_row(row)

    tab.sort(('target', 'filter', 'date'))
    tab.write('phot.txt', format='ascii.fixed_width_two_line',
              overwrite=True)

# now bin the data


def grouper(row):
    # bin every time the target, site, or filter changes and when Δt >= 1 hr
    return (row['target'], os.path.basename(row['file'])[:8], row['filter'], str(row['tbin']))


# resort for better grouping
tab.sort(('target', 'filter', 'file'))

tab['tbin'] = np.cumsum(
    np.r_[0, np.diff(Time(tab['date']).mjd * 24).astype(int)]
)
rows = []
for source, group in groupby(tab, grouper):
    g = vstack(list(group))
    i = (g['merr5'] < 0.2) * np.isfinite(g['m5'])

    # anything more than 3 pix away is probably bad
    if len(g[i]) > 2:
        # if there are more than 2 detections, use the median as a reference
        cm = np.median(g['dc'][i])
    else:
        # otherwise use the ephemeris position
        cm = 0
    i = i * (np.abs(g['dc'].data - cm) < 3)

    if not any(i):
        continue

    t = Time(g['date'][i])
    jd = t.jd.mean()
    fracday = t.mjd.mean() % 1
    tmtp = g['tmtp'].mean()
    rh = g['rh'][i].mean()
    delta = g['delta'][i].mean()
    phase = g['phase'][i].mean()
    n = i.sum()
    exptime = g['exptime'][i].sum()
    airmass = g['airmass'][i].mean()
    seeing = g['seeing'][i].mean()
    rho10k = g['rho10k'][i].mean()
    m = []
    merr = []
    for r in rho_labels:
        _m, sw = np.average(
            g['m{}'.format(r)][i],
            weights=g['merr{}'.format(r)][i]**-2,
            returned=True)
        merr.append(sw**-0.5)
        m.append(_m)

    rows.append([source[0], source[1], Time(jd, format='jd').iso[:-4],
                 fracday, tmtp, rh, delta, phase, source[2], g['catalog filter'][0],
                 n, exptime, airmass, seeing, rho10k] + m + merr)

del tab['tbin']


names = ['target', 'source', 'date', 'fracday', 'tmtp', 'rh', 'delta',
         'phase', 'filter', 'catalog filter', 'N exp', 'exptime',
         'airmass', 'seeing', 'rho10k']
names.extend(tab.colnames[-len(rho_labels) * 2:])
binned = Table(rows=rows, names=names, meta=tab.meta)


# color calculations and outburst search
colors = []
avg_colors = []
ostat = np.zeros(len(binned))
for target in set(binned['target']):
    i = binned['target'] == target
    # prev = binned[
    #     (binned['target'] == row['target'])
    #     #* (binned['filter'] == row['filter'])
    #     #* (binned['tmtp'] < row['tmtp'])
    #     #* (binned['tmtp'] > (row['tmtp'] - 14))
    # ]
    eph = Ephem.from_dict({
        'date': Time(binned['date'][i]),
        'rh': binned['rh'][i] * u.au,
        'delta': binned['delta'][i] * u.au,
        'phase': binned['phase'][i] * u.deg,
    })
    trend = CometaryTrends(eph, binned[f'm{color_rho}'][i] * u.mag,
                           binned[f'merr{color_rho}'][i] * u.mag,
                           filt=binned['catalog filter'][i],
                           color_transform=True,
                           logger=logger)

    gmr = trend.color('g', 'r')
    row = {
        'target': target
    }
    if gmr is not None:
        for j in range(len(gmr.c)):
            colors.append({
                'target': target,
                'date': gmr.t[j].iso,
                'filters': 'g-r',
                'c': gmr.c[j],
                'unc': gmr.c_unc[j]
            })

        row['n g-r'] = len(gmr.c)
        row['g-r'] = gmr.avg
        row['g-r unc'] = gmr.avg_unc

    rmi = trend.color('r', 'i')
    if rmi is not None:
        for j in range(len(rmi.c)):
            colors.append({
                'target': target,
                'date': rmi.t[j].iso,
                'filters': 'r-i',
                'c': rmi.c[j],
                'unc': rmi.c_unc[j]
            })

        row['n r-i'] = len(rmi.c)
        row['r-i'] = rmi.avg
        row['r-i unc'] = rmi.avg_unc

    if not (gmr is None and rmi is None):
        avg_colors.append(row)

    ostat[i] = trend.ostat(Phi=phase_HalleyMarcus, fixed_angular_size=True)

binned['ostat'] = ostat

colors = Table(colors)
avg_colors = Table(avg_colors)
for tab in (colors, avg_colors):
    tab.meta['comments'] = [f'colors measured with radius = {color_rho}']
colors.write('colors.txt', format='ascii.fixed_width_two_line',
             overwrite=True)
avg_colors.write('avg-colors.txt', format='ascii.fixed_width_two_line',
                 overwrite=True)

binned['fracday'].format = '{:.5f}'
binned['tmtp'].format = '{:.3f}'
binned['rh'].format = '{:.3f}'
binned['delta'].format = '{:.3f}'
binned['phase'].format = '{:.3f}'
binned['exptime'].format = '{:.0f}'
binned['airmass'].format = '{:.3f}'
binned['seeing'].format = '{:.2f}'
binned['rho10k'].format = '{:.2f}'
binned['ostat'].format = '{:.1f}'
for r in rho_labels:
    binned['m{}'.format(r)].format = '{:.3f}'
    binned['merr{}'.format(r)].format = '{:.3f}'

binned.write('phot-binned.txt', format='ascii.fixed_width_two_line',
             overwrite=True)
binned.sort('date', reverse=True)
binned.write('phot-binned-date-sort.txt', format='ascii.fixed_width_two_line',
             overwrite=True)

binned[binned['target'] == 'C/2014 UN271'].write(
    'c2014un271-look-phot-binned.txt',
    format='ascii.fixed_width_two_line',
    overwrite=True)

if len(new_objects) > 0:
    logger.info('New observations of %s', new_objects)
else:
    logger.info('No new observations')
