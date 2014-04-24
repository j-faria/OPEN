# In order to plot the GLS and BLS together 'run -i' this script:
magic = get_ipython().magic

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from OPEN.ext.get_rvN import get_rvn

import sys
sys.path.append('/home/joao/Work/Python')
from make_like_R import rstyle, pretify

set1 = pretify(1, set_rcParams=False)

colors = set1 # lets hope for less than 9 data-sets
t, rv, err = default.time, default.vrad, default.error # temporaries
tt = default.get_time_to_plot()
vel = np.zeros_like(tt)

## best solution found
best_ind = np.argmax(m.chains[1, :])
par_best = m.chains[2:, best_ind]

f = plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
ax = [plt.subplot(gs[0])]
ax.append(plt.subplot(gs[1], sharex=ax[0]))
# , ax = plt.subplots(nrows=2, sharex=True)


#### plot best solution
args = [tt] + list(par_best) + [vel]
get_rvn(*args)
ax[0].plot(tt, (vel - par_best[5])*1e3, '-', lw=1.5, color='k')
limits = ax[0].axis()

#### plot posterior predictive
# n_posterior_samples = 1000
# par_samples = sample(m.chains[2:, :].T, n_posterior_samples)
# for p in par_samples:
#     # if abs(p[4] - 2493817.60) > 1000: continue
#     args = [tt] + list(p) + [vel]
#     get_rvn(*args)
#     ax[0].plot(tt, (vel - p[5])*1e3, '-', lw=1.5, alpha=0.05)

# we have to do it again...
args = [tt] + list(par_best) + [vel]
get_rvn(*args)
ax[0].plot(tt, (vel - par_best[5])*1e3, '-', lw=1.5, color='k')



#### plot each files' values
for i, (fname, [n, nout]) in enumerate(sorted(default.provenance.iteritems())):
    mm = n-nout # how many values are there after restriction
    
    # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
    #                     height=err[:m], beam=0.5,\
    #                     pen=pg.mkPen(None))
                        # pen={'color': 0.8, 'width': 2})
    # p.addItem(e)
    # p.plot(t[:m], rv[:m], symbol='o')
    ax[0].errorbar(t[:mm], (rv[:mm]-par_best[5])*1e3, yerr=err[:mm]*1e3, \
                 fmt='o', color=colors[i], capsize=0)
    t, rv, err = t[mm:], rv[mm:], err[mm:]

ax[0].set_ylabel('RV [m/s]')
ax[0].set_ylim(list(limits[2:4]))


#### plot the residuals
t, rv, err = default.time, default.vrad, default.error # temporaries
vel = np.zeros_like(t)
args = [t] + list(par_best) + [vel]
get_rvn(*args)
for i, (fname, [n, nout]) in enumerate(sorted(default.provenance.iteritems())):
    mm = n-nout # how many values are there after restriction
    ax[1].errorbar(t[:mm], (rv[:mm]-vel[:mm])*1e3, yerr=err[:mm]*1e3, \
                 fmt='o', color=colors[i], capsize=0)
    t, rv, err, vel = t[mm:], rv[mm:], err[mm:], vel[:mm]
ax[1].axhline(y=0, lw=1.5, ls='--', color='k')

ax[1].set_xlabel('Time [days]')
ax[1].set_ylabel('O-C [m/s]')
plt.tight_layout()
plt.ticklabel_format(useOffset=False)

rstyle(ax, grid_toggle=True)
# ax[0].set_xticks([])
plt.show()