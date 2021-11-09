import warnings
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import ascii
from astropy.time import Time
from sbpy.data import Ephem

warnings.simplefilter('ignore', UserWarning)

tp = 2462889.37076

phot = ascii.read('phot.txt')
binned = ascii.read('phot-binned.txt')

phot = phot[(phot['target'] == 'C/2014 UN271')]
binned = binned[(binned['target'] == 'C/2014 UN271')]
bdate = Time(binned['date'])
date = Time(phot['date'])


alpha = 0.5
elinewidth = 1
aperture = '5'

fig = plt.figure(1, (10, 6))
fig.clear()

ax = fig.add_subplot(111)
ax.plot_date(date.plot_date, phot[f'm{aperture}'], ls='none', ms=2,
             label='LOOK singletons', color='tab:gray', marker='.')

ax.errorbar(bdate.plot_date, binned[f'm{aperture}'],
            binned[f'merr{aperture}'], ecolor='tab:gray', ls='none',
            elinewidth=elinewidth)

rax1 = ax.twinx()
eph = Ephem.from_horizons(
    'C/2014 UN271',
    epochs={
        'start': Time('2021-06-15'),
        'stop': Time('2021-08-01'),
        'step': 1 * u.d
    }
)
rax1.plot(eph['date'].plot_date, eph['rh'], label='rh')
rax1.plot(eph['date'].plot_date, eph['delta'], label='Δ')

rax2 = ax.twinx()
rax2.plot(eph['date'].plot_date, eph['phase'], label='Phase angle',
          color='tab:green')
rax2.plot(eph['date'].plot_date, eph['OrbPlaneAng'],
          label='Orbital plane angle',
          color='tab:red')


opts = {
    'r': {'color': 'tab:orange', 'marker': 's'},
    'g': {'color': 'tab:green', 'marker': 'o'},
    'i': {'color': 'tab:red', 'marker': '^'},
}
for filt in opts:
    i = binned['catalog filter'] == filt
    ax.plot_date(bdate.plot_date[i], binned[f'm{aperture}'][i], ls='none',
                 ms=4, label='LOOK ' + filt,
                 alpha=alpha, **opts[filt])

plt.setp(ax, ylabel='m (mag)', xlabel='Date (UTC)', ylim=[20.5, 18])
plt.setp(rax1, ylabel='rh / Δ (au)')
plt.setp(rax2, ylabel='Phase angle (deg)')
rax2.spines['right'].set_position(("axes", 1.15))

ax.minorticks_on()
rax1.minorticks_on()
rax2.minorticks_on()
fig.legend(loc='center right', frameon=False)
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
plt.tight_layout(rect=(0, 0, 0.825, 1))
plt.savefig('c2014un271-look.png')
