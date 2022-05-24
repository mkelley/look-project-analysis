import argparse
import numpy as np
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from sbpy.data import Ephem
from sbpy.activity import phase_HalleyMarcus
from mskpy.photometry import CometaryTrends
from lco_util import color_rho, setup_logger

parser = argparse.ArgumentParser()
parser.add_argument('-v', dest='level', action='store_const',
                    const='DEBUG', help='verbose mode')
parser.add_argument('-q', dest='level', action='store_const',
                    const='WARNING', help='quiet mode')
args = parser.parse_args()

logger = setup_logger(__file__, args.level)
logger.info('Start')

binned = ascii.read('phot-binned.txt')

# color calculations and outburst search
colors = []
avg_colors = []
ostat = np.zeros(len(binned))
for target in set(binned['target']):
    logger.debug(target)
    i = binned['target'] == target
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
                           logger=logger, log_level="DEBUG")

    logger.debug('g-r')
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

    logger.debug('r-i')
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

    ostat[i] = trend.ostat(
        k=4, dt=14, Phi=phase_HalleyMarcus, fixed_angular_size=True)

binned['ostat'] = ostat
binned.write('phot-binned.txt', format='ascii.fixed_width_two_line',
             overwrite=True)

colors = Table(colors)
for col in ['c', 'unc']:
    colors[col].format = '{:.3f}'

avg_colors = Table(avg_colors)
for col in ['g-r', 'g-r unc', 'r-i', 'r-i unc']:
    avg_colors[col].format = '{:.3f}'

for tab in (colors, avg_colors):
    tab.meta['comments'] = [f'colors measured with radius = {color_rho}']

colors.write('colors.txt', format='ascii.fixed_width_two_line',
             overwrite=True)
avg_colors.write('avg-colors.txt', format='ascii.fixed_width_two_line',
                 overwrite=True)

logger.info('Finish')
