import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.time import Time

# JPL K202/2, integrated to perihelion
# 2459198.500000000 = A.D. 2020-Dec-15 00:00:00.0000 TDB
#  EC= 7.356629428158953E-01 QR= 8.079228650349495E-01 IN= 1.394471268904226E+01
#  OM= 2.418511597598379E+02 W = 1.535440269626398E+02 Tp=  2459199.770713313017
#  N = 1.844532454935002E-01 MA= 3.597656128053511E+02 TA= 3.577283938231034E+02
#  A = 3.056411664870168E+00 AD= 5.304900464705386E+00 PR= 1.951714099889263E+03
tp = 2459199.770713313017

phot = ascii.read('phot.txt')
binned = ascii.read('phot-binned.txt')

target = 'C/2020 R4'
phot = phot[(phot['target'] == target)
            * (phot['filter'] == 'rp')]
binned = binned[(binned['target'] == target)
                * (binned['filter'] == 'rp')]

alpha = 0.5
elinewidth = 1

plt.clf()
ax = plt.gca()
ax.plot_date(Time(binned['date']).plot_date, binned['m5'], ls='none',
             ms=4, label='LOOK', marker='s', color='tab:pink', alpha=alpha)
ax.errorbar(Time(binned['date']).plot_date, binned['m5'],
            binned['merr5'], ecolor='tab:gray', ls='none', elinewidth=elinewidth)

ax.plot_date(Time(phot['date']).plot_date, phot['m5'], ls='none', ms=2,
             label='LOOK singletons', color='tab:gray', marker='.')

plt.setp(ax, ylabel='m (mag)', xlabel='Date (UTC)')
ax.minorticks_on()
ax.invert_yaxis()
ax.legend(loc='lower left')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f'{target.replace(" ", "").replace("/", "").lower()}-look.png')
