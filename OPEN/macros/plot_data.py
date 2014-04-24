# 'run -i' this script:
magic = get_ipython().magic

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import patches
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


f = plt.figure(dpi=150, figsize=(9,5))
# gs = gridspec.GridSpec(2, 1)
# ax = [plt.subplot(gs[0]), plt.subplot(gs[1])]
# plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.15, hspace=0.4)
ax = plt.subplot(111)
plt.subplots_adjust(bottom=0.15)

#### plot each files' values
for i, (fname, [n, nout]) in enumerate(sorted(default.provenance.iteritems())):
    mm = n-nout # how many values are there after restriction
    ax.errorbar(t[:mm], rv[:mm]*1e3, yerr=err[:mm]*1e3, \
                 fmt='o', color=colors[i], capsize=0, alpha=0.8)
    t, rv, err = t[mm:], rv[mm:], err[mm:]

ax.set_xlabel('BJD - 2400000 [days]')
ax.set_ylabel('RV [m/s]')
ax.set_title('HD41248', loc='right')
plt.ticklabel_format(useOffset=False)
# ax[0].set_xlim([49000, 54000])


do_rectangle = True
if do_rectangle:
	rect1 = patches.Rectangle((52800,-13), 3000, 26, color=colors[1], alpha=0.4)
	ax.add_patch(rect1)
	ax.text(54500, 11.5, 'Jenkins et al. (2013)')


rstyle([ax], grid_toggle=True)

plt.show()