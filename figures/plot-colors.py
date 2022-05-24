#from glob import glob
import numpy as np
#import scipy.ndimage as nd
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
#import astropy.units as u
#from astropy.io import fits
from astropy.io import ascii
#from astropy.table import Table
#from astropy.time import Time
#from astropy.coordinates import SkyCoord
from astropy.visualization import hist
#import mskpy

colors = ascii.read('colors.txt')
gmr = colors['filters'] == 'g-r'
targets = set(colors['target'])
weights = colors['unc']**-2

avg_colors = ascii.read('avg-colors.txt')

fig = plt.figure(1, clear=True)
ax = fig.gca()

n, bins, patches = hist(colors['c'][gmr], bins='freedman',
                        color='0.75')

x = np.linspace(0.2, 0.9, 300)
top = ax.get_ylim()[1]
for i, target in enumerate(targets):
    # hist(colors[gmr * (colors['target'] == target)])
    j = gmr * (colors['target'] == target)
    if j.sum() > 8:
        kernel = gaussian_kde(colors['c'][j], weights=weights[j])
        line, = ax.plot(x, kernel(x))

        c = avg_colors[avg_colors['target'] == target][0]
        ax.axvline(c['g-r'], color=line.get_color(), lw=1)
        ax.text(c['g-r'] - 0.003,  # top * ((i * 0.1 % 0.8) + 0.02),
                kernel(x).max() + 0.02 * top,
                target, rotation=90, ha='center', va='bottom',
                size=8)


ax.minorticks_on()
plt.setp(ax, ylabel='Number', xlabel='g-r (mag)', xlim=(0.2, 0.9))
plt.tight_layout(pad=0.4)
fig.savefig('color-histogram.png', dpi=200)
