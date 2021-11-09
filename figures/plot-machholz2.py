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

phot = phot[(phot['target'] == '141P')
            * (phot['filter'] == 'rp')]
binned = binned[(binned['target'] == '141P')
                * (binned['filter'] == 'rp')]

ztf = ascii.read('/home/msk/data/targets/141p/analysis/outburst-mar21/'
                 'rap5/ztf-phot-rap5-final-binned.txt')
ztf = ztf[ztf['date'] > '2021-02-15']
#ztf[ztf['filter'] == 'ZTF_g']['m'] -= 0.28

ldt = ascii.read('''
date,rap,m,m err
#2021-03-03 05:05,12,14.13,0.02
2021-03-03 05:05,5,15.05,0.02
2021-03-05 02:40,5,16.72,0.02
''', delimiter=',')

alpha = 0.5
elinewidth = 1

plt.clf()
ax = plt.gca()
for filt, color in zip(('ZTF_r', 'ZTF_g'), (0, 0.28)):
    i = ztf['filter'] == filt
    label = f'ZTF {filt[-1]} '
    if color > 0:
        label += f'{color:+.2f} mag'

    ax.plot_date(Time(ztf['date']).plot_date[i], ztf['m'][i] - color, ls='none',
                 ms=4, label=label, alpha=alpha)
    ax.errorbar(Time(ztf['date']).plot_date[i], ztf['m'][i] - color,
                ztf['m err'][i], ecolor='tab:gray', ls='none', elinewidth=elinewidth)

ax.plot_date(Time(ldt['date']).plot_date, ldt['m'], ls='none',
             ms=4, label='LDT', color='tab:purple', alpha=alpha)
ax.errorbar(Time(ldt['date']).plot_date, ldt['m'],
            ldt['m err'], ecolor='tab:gray', ls='none', elinewidth=elinewidth)

ax.plot_date(Time(binned['date']).plot_date, binned['m5'], ls='none',
             ms=4, label='LOOK', marker='s', color='tab:pink', alpha=alpha)
ax.errorbar(Time(binned['date']).plot_date, binned['m5'],
            binned['merr5'], ecolor='tab:gray', ls='none', elinewidth=elinewidth)

ax.plot_date(Time(phot['date']).plot_date, phot['m5'], ls='none', ms=2,
             label='LOOK singletons', color='tab:gray', marker='.')

plt.setp(ax, ylabel='m (mag)', xlabel='Date (UTC)')
ax.minorticks_on()
ax.invert_yaxis()
ax.legend(loc='upper left')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
plt.tight_layout()
plt.savefig('141p-look.png')
