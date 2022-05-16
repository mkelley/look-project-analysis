import os
import json
import logging
import argparse
from glob import glob
from collections import defaultdict
from itertools import groupby
import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
from astropy.io import ascii, fits
from astropy.time import Time
from astropy.table import Table
from astropy.visualization import imshow_norm, ZScaleInterval, LogStretch
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from ccdproc import Combiner
from sbpy.activity import phase_HalleyMarcus as Phi
import lco_util

parser = argparse.ArgumentParser()
parser.add_argument('--window', default=14, help='use this time window (days)')
parser.add_argument('--target', help='only plot this target')
parser.add_argument('--stack', help='replot this stack file (implies --force)')
parser.add_argument('--force', '-f', action='store_true',
                    help='overwrite existing images')
parser.add_argument('-v', action='store_true', help='verbose mode')
parser.add_argument('-q', action='store_true', help='quiet mode')
args = parser.parse_args()

logger = logging.getLogger(__name__)
if len(logger.handlers) == 0:
    logger.addHandler(logging.StreamHandler())
    logger.addHandler(logging.FileHandler('look-project-analysis.log'))

if args.v:
    logger.setLevel('DEBUG')
elif args.q:
    logger.setLevel('WARNING')
else:
    logger.setLevel('INFO')

# invert file -> cluster associations, but only track "median" files.
with open('stack-clusters.json', 'r') as inf:
    associations = json.load(inf)
stacks = defaultdict(set)
for fn, stack_files in associations.items():
    for stack in stack_files:
        if stack.endswith('med.fits'):
            stacks[stack].add(fn)


def grouper(stack):
    return fits.getheader(stack[0])['target']

half_box = 230
y0, x0 = 537, 537
cutout = np.s_[y0 - half_box: y0 + half_box + 1,
               x0 - half_box: x0 + half_box + 1]

for target, group in groupby(sorted(stacks.items()), key=grouper):
    if args.target is not None:
        if target != args.target:
            continue

    group = list(group)
    logger.info('%s %d stacks', target, len(group))

    headers = [fits.getheader(stack[0]) for stack in group]
    obsdates = np.array([h['midmjd'] for h in headers])
    filters = np.array([h['filter'] for h in headers])

    for i, (stack_file, sources) in enumerate(group):        
        plot_file = stack_file.replace('.fits', '.png')
        logger.debug('checking %s', plot_file)
        if os.path.exists(plot_file) and not args.force:
            # if stack file creation time is older than plot, then skip this stack
            if os.stat(stack_file).st_ctime < os.stat(plot_file).st_ctime:
                continue

        # else this plot needs to be updated
        logger.info('updating %s', plot_file)

        nightly = fits.getdata(stack_file)
        h0 = headers[i]

        # find recent observations with the same filter, if any
        dt = obsdates - obsdates[i]
        last_two_weeks = (dt < 0) * (dt > -args.window) * (filters == filters[i])
        if any(last_two_weeks):
            two_weeks = []
            for f in [group[j][0] for j in np.flatnonzero(last_two_weeks)]:
                im, h = fits.getdata(f, header=True)
                if h['mzpunc'] > 0.1:
                    continue

                two_weeks.append(
                    im * (h['rh'] / h0['rh']) ** 4 * h['delta'] / h0['delta']
                    * (Phi(h['phase']) / Phi(h0['phase']))**-1
                )

            two_weeks = CCDData(two_weeks, unit='adu')
            combiner = Combiner(two_weeks)

            med = combiner.median_combine().data
            diff = nightly - med
        else:
            med = np.zeros_like(nightly)
            diff = med

        fig = plt.figure(1, (8, 16 / 3), clear=True)
        tax = fig.add_axes((0, 0.5, 1, 0.5))
        bax = fig.add_axes((0, 0, 1, 0.5))

        mean, median, stdev = sigma_clipped_stats(nightly)
        opts = {
            'vmin': median,
            'vmax': np.mean(nightly[
                int(h0['crpix2']) - 6:int(h0['crpix2']) + 5,
                int(h0['crpix1']) - 6:int(h0['crpix1']) + 5
            ])
        }

        im = np.c_[nightly[cutout], med[cutout], diff[cutout]]
        imshow_norm(im, tax, interval=ZScaleInterval())
        imshow_norm(im, bax, stretch=LogStretch(), **opts)

        x = np.array((0.4, 0.45, 0.55, 0.6)) * im.shape[0]
        opts = {'color': '0.75', 'lw': 1}
        for ax in (tax, bax):
            for dx in np.r_[0, 1, 2] * im.shape[0]:
                ax.plot(x[:2] + dx, x[:2], **opts)
                ax.plot(x[:2][::-1] + dx, x[2:], **opts)
                ax.plot(x[2:] + dx, x[2:], **opts)
                ax.plot(x[2:][::-1] + dx, x[:2], **opts)

        plt.setp((tax, bax), frame_on=False, xticks=[], yticks=[])
        plt.savefig(plot_file, dpi=200)
