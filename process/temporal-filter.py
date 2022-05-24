import os
import json
import logging
import argparse
from collections import defaultdict
from itertools import groupby
import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.visualization import imshow_norm, ZScaleInterval, LogStretch
from astropy.stats import sigma_clipped_stats
from ccdproc import Combiner
from sbpy.activity import phase_HalleyMarcus as Phi
import lco_util

parser = argparse.ArgumentParser()
parser.add_argument('--window', default=14, help='use this time window (days)')
parser.add_argument('--target', help='only plot this target')
parser.add_argument('--stack', help='replot this stack file (implies --force)')
parser.add_argument('--force', '-f', action='store_true',
                    help='overwrite existing images')
parser.add_argument('-v', dest='level', action='store_const',
                    const='DEBUG', help='console verbose mode')
parser.add_argument('-q', dest='level', action='store_const',
                    const='WARNING', help='console quiet mode')
args = parser.parse_args()

logger = lco_util.setup_logger(__file__, args.level)
logger.info('Start')

formatter = logging.Formatter(
    '%(name)s : %(levelname)s : %(asctime)s : %(message)s')
for h in logger.handlers:
    h.setFormatter(formatter)

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
    logger.debug('%s %d stack%s', target, len(
        group), '' if len(group) == 1 else 's')

    headers = [fits.getheader(stack[0]) for stack in group]
    obsdates = np.array([h['midmjd'] for h in headers])
    filters = np.array([h['filter'] for h in headers])

    for i, (stack_file, sources) in enumerate(group):
        filtered_file = (stack_file
                         .replace('stacks/', 'processed/')
                         .replace('.fits', '_tf.fits'))
        plot_file = (filtered_file
                     .replace('.fits', '.jpg')
                     .replace('processed/', 'plots/'))

        for f in (filtered_file, plot_file):
            dirname = os.path.dirname(f)
            if not os.path.exists(dirname):
                os.system(f'mkdir {dirname}')

        logger.debug('checking %s', plot_file)
        if os.path.exists(filtered_file) and not args.force:
            # if stack file creation time is older than filtered data, then skip
            # this stack
            if os.stat(stack_file).st_ctime < os.stat(filtered_file).st_ctime:
                continue

        # else this plot needs to be updated
        logger.info('updating %s', filtered_file)

        nightly = fits.getdata(stack_file)
        h0 = headers[i]
        h0["window"] = (args.window, 'time window for baseline files (days)')

        # find recent observations with the same filter, if any
        dt = obsdates - obsdates[i]
        last_two_weeks = (dt < 0) * (dt > -args.window) * \
            (filters == filters[i])
        if any(last_two_weeks):
            two_weeks = []
            files = []
            for f in [group[j][0] for j in np.flatnonzero(last_two_weeks)]:
                files.append(f)
                im, h = fits.getdata(f, header=True)

                two_weeks.append(
                    im * (h['rh'] / h0['rh']) ** 4 * h['delta'] / h0['delta']
                    * (Phi(h['phase']) / Phi(h0['phase']))**-1
                )

            h0["nbl"] = (len(files), 'number of baseline files')
            h0.add_comment(f"Baseline files: {str(files)}")

            two_weeks = CCDData(two_weeks, unit='adu')
            combiner = Combiner(two_weeks)

            med = combiner.median_combine().data
            diff = nightly - med
        else:
            med = np.zeros_like(nightly)
            diff = med

        hdu = fits.HDUList()
        hdu.append(fits.ImageHDU(diff, h0, name="filtered"))
        hdu.append(fits.ImageHDU(nightly, name="nightly"))
        hdu.append(fits.ImageHDU(med, name="baseline"))
        hdu.writeto(filtered_file, overwrite=True)

        fig = plt.figure(1, (6, 12 / 3), clear=True)
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
        plt.savefig(plot_file, dpi=100)

logger.info('Finish')
