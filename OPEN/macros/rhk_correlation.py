# -*- coding: utf-8 -*-
# 'run -i' this script:
magic = get_ipython().magic

from OPEN.shell_colors import blue
from OPEN.classes import params
from scipy.odr import odrpack
import sys
import matplotlib.ticker as ticker
import re
import os

try:
	save = sys.argv[1]
except IndexError:
	save = None

# this is Seaborn's "colorblind" pallete
colors = ['#0072b2', '#009e73', '#d55e00', '#cc79a7', '#f0e442', '#56b4e9']
# this is Seaborn's "muted" pallete
# colors = ['#4878cf', '#6acc65', '#d65f5f', '#b47cc7', '#c4ad66', '#77bedb']


full_path = default.provenance.keys()[0]
bn = os.path.basename(full_path)
i = bn.rfind('_harps_mean_corr.rdb')
if i == -1:
    i = bn.rfind('_harps_mean.rdb')
star = bn[:i]



pr = sp.stats.pearsonr(default.extras.rhk, default.vrad)
sr = sp.stats.spearmanr(default.extras.rhk, default.vrad)

print blue('[Pearson correlation]') + ' r=%f, p-value=%f' % pr
print blue('[Spearman correlation]') + ' r=%f, p-value=%f' % sr

v1 = default.extras.rhk
e1 = default.extras.sig_rhk
v2 = default.vrad
e2 = default.error

with plt.rc_context(params):

    class MyFormatter(ticker.ScalarFormatter):
        def __call__(self, x, pos=None):
            # call the original LogFormatter
            rv = ticker.ScalarFormatter.__call__(self, x, pos)
            # check if we really use TeX
            if plt.rcParams["text.usetex"]:
                # if we have the string ^{- there is a negative exponent
                # where the minus sign is replaced by the short hyphen
                rv = re.sub(r'-', r'\mhyphen', rv)
            return rv

    figwidth = 3.543311946  # in inches = \hsize = 256.0748pt
    figheight = 0.75 * figwidth

    fig = plt.figure(figsize=(figwidth, figheight))
    ax = fig.add_subplot(111)
    ax.set_title(star, loc='right', fontsize=params['font.size'])
    
    (_, caps, barlines) = ax.errorbar(v1, v2, xerr=e1, yerr=e2,
    	                              fmt='o', color='k', mec='none', ms=2, capsize=0, elinewidth=0.5)

    for cap in caps:
        cap.set_alpha(0.3)
    for barline in barlines:
        barline.set_alpha(0.3)

    # non-weighted fit, OLS
    m, b = np.polyfit(v1, v2, 1)
    yp = np.polyval([m, b], v1)
    # ax.plot(v1, yp, '-k')

    # weghted fit, ODR (only do if errors on variables)
    if not (np.count_nonzero(e1) == 0 or np.count_nonzero(e2) == 0):
        def f(B, x): 
        	return B[0]*x + B[1]
        linear = odrpack.Model(f)
        data = odrpack.Data(v1, v2, wd=1./e1**2, we=1./e2**2)
        odr = odrpack.ODR(data, linear, beta0=[1., 1.])
        output = odr.run()
        yp = np.polyval(output.beta, v1)
        # ax.plot(v1, yp, '--k')

    ax.set_ylabel(r'RV [$\ms$]')
    ax.set_xlabel(r"$\log R'_{HK}$")
    ax.xaxis.set_major_formatter(MyFormatter())
    fig.tight_layout()

    if save is not None:
    	fig.savefig(save)

    plt.show()