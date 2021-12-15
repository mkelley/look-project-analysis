import warnings
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import ascii
from astropy.time import Time
from sbpy.data import Ephem
from lco_util import target_to_filename

warnings.simplefilter('ignore', UserWarning)

tp = 2459520.55679
target = '67P'

phot = ascii.read('phot.txt')
binned = ascii.read('phot-binned.txt')

phot = phot[(phot['target'] == target)]
binned = binned[(binned['target'] == target)]
bdate = Time(binned['date'])
date = Time(phot['date'])

alpha = 0.5
elinewidth = 1
aperture = '5'
color = {
    'g': 0.62,
    'r': 0,
    'i': -0.26
}

fig = plt.figure(1, (10, 6))
fig.clear()

ax = fig.add_subplot(111)
c = np.array([color[filt] for filt in phot['catalog filter']])
ax.plot_date(date.plot_date, phot[f'm{aperture}'] - c, ls='none', ms=2,
             label='LOOK singletons', color='tab:gray', marker='.')

opts = {
    'r': {'color': 'tab:orange', 'marker': 's'},
    'g': {'color': 'tab:green', 'marker': 'o'},
    'i': {'color': 'tab:red', 'marker': '^'},
}
for filt in opts:
    i = (binned['catalog filter'] == filt) * (binned['date'] > '2021-10-01')
    if sum(i) == 0:
        continue
    ax.plot_date(bdate.plot_date[i], binned[f'm{aperture}'][i] - color[filt],
                 ls='none', ms=4, label='LOOK ' + filt,
                 alpha=alpha, **opts[filt])

    ax.errorbar(bdate.plot_date[i], binned[f'm{aperture}'][i] - color[filt],
                binned[f'merr{aperture}'][i], ecolor='tab:gray', ls='none',
                elinewidth=elinewidth)

plt.setp(ax, ylabel=f'$m_{{{aperture}}}$ (mag)', xlabel='Date (UTC)', ylim=[13, 11],
         title=target)

ax.minorticks_on()
fig.legend(loc='center right', frameon=False)
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
plt.tight_layout(rect=(0, 0, 0.825, 1))
plt.savefig(f'{target_to_filename(target)}-look.png')
