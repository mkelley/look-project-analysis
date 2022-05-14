#!/usr/bin/env python3
import os
import logging
import argparse
import warnings
from itertools import groupby
from astropy.nddata import CCDData
import numpy as np
from scipy.cluster import hierarchy
from astropy.time import Time
from astropy.io import ascii, fits
from astropy.table import vstack
from astropy.wcs import WCS, FITSFixedWarning
import astropy.units as u
from reproject import reproject_exact, reproject_adaptive, reproject_interp
from ccdproc import Combiner
from sbpy.data import Ephem, Names, TimeScaleWarning
from sbpy.data.names import TargetNameParseError
#from mskpy.image import align_by_offset
from lco_util import (assumed_gmr, locations,
                      filters as filter_to_ps1, target_to_filename)

parser = argparse.ArgumentParser()
parser.add_argument('--target', help='only stack this target')
parser.add_argument('--date', help='only stack on this UT date')
parser.add_argument('--cluster-dt', type=float, default=1,
                    help='distance threshold for clustering (hours)')
parser.add_argument('--filter', help='only stack this filter')
parser.add_argument('--force', '-f', action='store_true', help='overwrite existing images')
parser.add_argument('--size', '-s', default=1000, type=int, help='image dimensions')
parser.add_argument('-v', action='store_true', help='verbose mode')
parser.add_argument('-q', action='store_true', help='quiet mode')
args = parser.parse_args()

warnings.simplefilter(
    'ignore', (UserWarning, FITSFixedWarning, TimeScaleWarning))

level = (2 + int(args.q) - int(args.v)) * 10
logging.basicConfig(level=level)

phot = ascii.read('phot.txt')
phot['mjd'] = Time(phot['date']).mjd
phot['site'] = [f.split('/')[2].split('-')[0] for f in phot['file']]

phot.sort('date')
phot.sort(('target', 'site', 'filter'))
targets = set(phot['target']) if args.target is None else [args.target]
filters = set(phot['filter']) if args.filter is None else [args.filter]

shape = np.array((args.size, args.size))


def grouper(row):
    return (row['target'], row['site'], row['filter'])


for (target, site, filter), group in groupby(phot, grouper):
    if (target not in targets) or (filter not in filters):
        continue

    group = vstack(list(group))

    target_fn = target_to_filename(target)
    os.system(f'mkdir -p stacks/{target_fn}')

    if len(group) > 1:
        clusters = hierarchy.fclusterdata(
            group['mjd'].data[:, np.newaxis],
            args.cluster_dt / 24,
            criterion='distance'
        )
    else:
        clusters = np.array([1])

    for cluster in [group[clusters == c] for c in set(clusters)]:
        cluster.sort('date')
        dates = Time(cluster['mjd'], format='mjd')
        date = Time(dates.mjd.mean(), format='mjd')

        if args.date is not None and date.iso[:10] != args.date:
            continue

        date_fn = date.isot.replace('-', '').replace(':', '')[:13]
        prefix = f"stacks/{target_fn}/{target_fn}_{date_fn}_{cluster['rh'].mean():.3f}au_{filter}_look"

        if all([os.path.exists(f'{prefix}_{label}.fits') for label in ('avg', 'med')]) and not args.force:
            continue

        logging.debug(prefix)

        # create WCS objects in the comet's rest frame for each image and the stacked image
        opts = dict(epochs=dates, location=locations[site[:3]])
        try:
            target_type = Names.asteroid_or_comet(target)
        except TargetNameParseError:
            target_type = 'asteroid'

        if target_type == 'comet':
            id_type = 'designation'
            opts['closest_apparition'] = True
            opts['no_fragments'] = True
        else:
            id_type = 'smallbody'

        eph = Ephem.from_horizons(target, id_type, **opts)

        wcs = []
        for i, f in enumerate(cluster['file']):
            w = WCS(fits.getheader(f, 1))
            w.wcs.crval[0] += (eph['ra'][0] - eph['ra'][i]).to_value('deg')
            w.wcs.crval[1] += (eph['dec'][0] - eph['dec'][i]).to_value('deg')
            wcs.append(w)

        wcs0 = WCS()
        wcs0.wcs.radesys = 'ICRS'
        wcs0.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs0.wcs.crval = cluster['ra'][0], cluster['dec'][0]
        wcs0.wcs.crpix = shape // 2
        #ps = np.abs(wcs[0].wcs.cd).mean() * 2
        # ps0 = np.mean((
        #     np.hypot(*wcs[0].wcs.cd[0]),
        #     np.hypot(*wcs[0].wcs.cd[1]),
        # )) * 0.98
        # ps0 = u.Quantity(wcs[0].proj_plane_pixel_scales()).mean().value
        #wcs0.wcs.cdelt = -ps0, ps0
        wcs0.wcs.cd = -wcs[0].wcs.cd
        wcs0.wcs.dateobs = date.isot
        wcs0.wcs.mjdobs = date.mjd
        Omega0 = wcs0.proj_plane_pixel_area()

        # background subtract, photometric calibration, re-projection
        stack = np.ma.empty((len(cluster['file']), shape[0], shape[1]))
        stack.mask = np.zeros(stack.shape, bool)
        scale = np.empty(len(stack))
        for i in range(len(cluster['file'])):
            im0 = fits.getdata(cluster['file'][i])
            im0 -= fits.getdata(
                f"backgrounds/{os.path.basename(cluster['file'][i][:-3])}")
            mzp = cluster['zp'][i] + \
                cluster['color cor'][i] * assumed_gmr[target]
            scale[i] = 10**(-0.4 * (mzp - 25))

            # conserve flux
            Omega = wcs[i].proj_plane_pixel_area()
            scale[i] *= (Omega / Omega0).value

            im0 *= scale[i]

            # I'm not sure how best to avoid interpolation artifacts.  It may be because
            # the rotation is usually very small.  Using nearest-neighbor here looks good.
            # but bicubic is best for actual interpolations.
            # stack.data[i], cov = reproject_adaptive((im0, wcs[i]), wcs0, shape,
            #                                     order='nearest-neighbor')
            # ...
            # Currently not aligning with north and using "exact" (conserving flux)
            stack.data[i], cov = reproject_exact(
                (im0, wcs[i]), wcs0, shape, parallel=True)
            stack.mask[i] = cov == 0

        stack = CCDData(stack, unit='adu')
        combiner = Combiner(stack)

        avg = combiner.average_combine().data
        if len(stack.data) > 2:
            med = combiner.median_combine().data
        else:
            med = None

        hdu = wcs0.to_fits()
        hdu[0].header['target'] = target
        hdu[0].header['site'] = site
        hdu[0].header['filter'] = filter, 'observation filter'
        hdu[0].header['calfilt'] = filter_to_ps1[filter], 'data is calibrated to this PS1 filter'
        hdu[0].header['gain'] = scale.mean()**-1, 'mean eff. gain per image'
        hdu[0].header['ra'] = cluster['ra'].mean(), 'mean observed RA, deg'
        hdu[0].header['dec'] = cluster['dec'].mean(), 'mean observed Dec, deg'
        hdu[0].header['rh'] = cluster['rh'].mean()
        hdu[0].header['delta'] = cluster['delta'].mean()
        hdu[0].header['phase'] = cluster['phase'].mean()
        hdu[0].header['rho10k'] = cluster['rho10k'].mean()
        hdu[0].header['airmass'] = cluster['airmass'].mean()
        hdu[0].header['seeing'] = (cluster['seeing'].mean(),
                                   'mean approximate seeing, arcsec')

        try:
            hdu[0].header['tmtp'] = cluster['tmtp'].mean(), 'T-Tp'
        except ValueError:
            # in case tmtp is masked
            pass
        hdu[0].header['exptime'] = cluster['exptime'].sum()
        hdu[0].header['mzp'] = 25.0, 'magnitude zero point'
        hdu[0].header['mzpunc'] = (np.sqrt(np.sum(cluster['zp err']**2)) / np.sqrt(len(cluster)),
                                   'zeropoint uncertainty')
        hdu[0].header['colorcor'] = (cluster['color cor'].mean(),
                                     'mean color correction')
        hdu[0].header['tgtcolor'] = assumed_gmr[target], 'target g-r used for calibration'

        hdu[0].header['nimages'] = len(cluster), 'number of images combined'
        hdu[0].header.add_comment(
            'Files combined: {}'.format([os.path.basename(f)
                                         for f in cluster['file']]))

        for im, label in zip((avg, med), ('avg', 'med')):
            if im is None:
                continue
            hdu[0].header['combfunc'] = label, 'combination function'
            hdu[0].data = im.astype(np.float32)
            hdu.writeto(f'{prefix}_{label}.fits', overwrite=True)
