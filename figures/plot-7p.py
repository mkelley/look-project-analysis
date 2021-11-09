import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.time import Time
from astropy.timeseries import LombScargle

# JPL K202/2, integrated to perihelion
# 2459198.500000000 = A.D. 2020-Dec-15 00:00:00.0000 TDB
#  EC= 7.356629428158953E-01 QR= 8.079228650349495E-01 IN= 1.394471268904226E+01
#  OM= 2.418511597598379E+02 W = 1.535440269626398E+02 Tp=  2459199.770713313017
#  N = 1.844532454935002E-01 MA= 3.597656128053511E+02 TA= 3.577283938231034E+02
#  A = 3.056411664870168E+00 AD= 5.304900464705386E+00 PR= 1.951714099889263E+03
tp = 2459199.770713313017

phot = ascii.read('phot.txt')
binned = ascii.read('phot-binned.txt')

phot = phot[(phot['target'] == '7P')
            * (phot['filter'] == 'rp')
            * (phot['date'] > '2021-05-15')]
binned = binned[(binned['target'] == '7P')
                * (binned['filter'] == 'rp')
                * (binned['date'] > '2021-05-15')]
bdate = Time(binned['date'])
date = Time(phot['date'])

# On 2021 June 6.842 UTC, the comet had an apparent brightness of r = 13.48 +/- 0.06 mag. By our next epoch on 2021 June 7.858 UTC, the comet had brightened to r = 13.19 +/- 0.03 mag, suggesting an outburst strength of -0.29 +/- 0.07 mag. During the second epoch, the comet was at 1.24 AU from the Sun, and 0.44 AU from the Earth. The event was confirmed in Zwicky Transient Facility (ZTF; Bellm et al. 2019, PASP, 131, a8002) data, where the comet had an apparent brightness g = 13.68 +/- 0.03 mag on 2021 June 8.439 UTC. The g-r color of the comet is 0.34 mag, as measured from ZTF data on 2021 June 4.
atel14687 = 59372.858

alpha = 0.5
elinewidth = 1

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
ax.plot_date(bdate.plot_date, binned['m5'], ls='none',
             ms=4, label='LOOK', marker='s', color='tab:pink', alpha=alpha)
ax.errorbar(bdate.plot_date, binned['m5'],
            binned['merr5'], ecolor='tab:gray', ls='none', elinewidth=elinewidth)

ax.plot_date(date.plot_date, phot['m5'], ls='none', ms=2,
             label='LOOK singletons', color='tab:gray', marker='.')

plt.setp(ax, ylabel='m (mag)', xlabel='Date (UTC)',
         ylim=[16, min(ax.get_ylim())])
ax.minorticks_on()
ax.legend(loc='lower left')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
plt.tight_layout()
plt.savefig('7p-look.png')

i = binned['tmtp'] > 10
H = binned['m5'] - 5 * np.log10(binned['rh']) - 2.5 * np.log10(binned['delta'])

fig3 = plt.figure(3)
fig3.clf()
ax3 = fig3.add_subplot(111)
frequency, power = LombScargle(bdate.mjd[i], H[i]).autopower()
ax3.plot(frequency, power)
fig3.savefig('7p-look-periodigram.png')

P = 24 / 1.22

fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
ax2.errorbar((bdate.mjd[i] / (P / 24)) % 1, H[i],
             binned['merr5'][i], ecolor='tab:gray', ls='none',
             elinewidth=elinewidth, ms=4, marker='s',
             color='tab:pink', alpha=alpha, label='LOOK')
ax2.axvline(atel14687 / (P / 24) % 1, label='ATel #14687 (GIT)')
plt.setp(ax2, ylabel='H (mag)', xlabel=f'Phase (P={P} hr)')
ax2.minorticks_on()
ax2.invert_yaxis()
fig2.tight_layout()
fig2.savefig('7p-look-phased.png')
