# 'run -i' this script:
magic = get_ipython().magic

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

import sys
sys.path.append('/home/joao/Work/Python')
from make_like_R import rstyle, pretify

# magic('read 14her.rv -d')
# try: 
# 	default.per
# except:
# 	magic('per obs --ofac=10')

set1 = pretify(1, set_rcParams=False)

colors = set1 
t, rv, err = default.time, default.vrad, default.error # temporaries


f = plt.figure(dpi=120)
gs = gridspec.GridSpec(2, 1)
ax = [plt.subplot(gs[0]), plt.subplot(gs[1])]
plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.15, hspace=0.4)

#### plot each files' values
for i, (fname, [n, nout]) in enumerate(sorted(default.provenance.iteritems())):
    mm = n-nout # how many values are there after restriction
    ax[0].errorbar(t[:mm]-2400000, rv[:mm], yerr=err[:mm], \
                 fmt='o', color=colors[i], capsize=0)
    t, rv, err = t[mm:], rv[mm:], err[mm:]

ax[0].set_xlabel('BJD - 2400000 [days]')
ax[0].set_ylabel('RV [m/s]')
ax[0].set_title('14 Her', loc='right')
# ax[0].set_xlim([49000, 54000])



ax[1].semilogx(1/default.per.freq, default.per.power, 'k', lw=1.5)

ax[1].set_xlabel('Period [days]')
if default.per.name == 'Bayesian Lomb-Scargle':
	ax[1].set_ylabel(r'P($f$|D,I)')
else:
	ax[1].set_ylabel('Normalized Power')



rstyle(ax, grid_toggle=True)

plt.show()