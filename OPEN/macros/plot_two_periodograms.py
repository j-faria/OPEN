# In order to plot the periodogram of the RV and the FWHM together 'run -i' this script:
magic = get_ipython().magic

import sys
import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, IndexFormatter
# plt.rcParams['text.latex.preamble'].append(r'\usepackage{lmodern}')
# plt.rcParams["text.latex.preamble"].append(r'\mathchardef\mhyphen="2D')
params = {'text.latex.preamble': [r'\usepackage{lmodern}', r'\mathchardef\mhyphen="2D'],
          'text.usetex' : True,
          'font.size'   : 8,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          'axes.unicode_minus': True,
          }

# magic('per obs --bayes')  # to calculate GLS for RV
# magic('per fwhm --bayes')  # to calculate GLS for RV

magic('per obs --noplot')  # to calculate GLS for RV
per_rv = default.per
# magic('per resid --noplot --bfap')  # to calculate GLS for RV
per_resid = default.resid_per


figwidth = 3.543311946  # in inches = \hsize = 256.0748pt
figheight = 0.75 * figwidth


with plt.rc_context(params):
	fig = plt.figure(figsize=(figwidth, figheight))
	gs = plt.GridSpec(2, 1)
	ax = [plt.subplot(gs[0]), plt.subplot(gs[1])]
	# plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.15, hspace=0.4)

	# if per_rv.name == 'Bayesian Lomb-Scargle':
	# 	ax[0].set_ylabel(r'Normalized  P($f$ |D,I)')
	# 	m = max(per_rv.power)
	# 	ax[0].set_ylim([0,1.05])
	# else:
	ax[0].set_ylabel('Power')

	ax[0].semilogx(1./per_rv.freq, per_rv.power, 'k')
	ax[0].axhline(per_rv.powerFAP_1, c='k', lw=1, ls='--', dashes=(4,3), alpha=0.6, label=r'1\% FAP')
	ax[0].set_xlim([0.5, 20])
	ax[0].set_ylim([0, 20])
	plt.setp(ax[0].get_xticklabels(), visible=False)
	ax[0].legend(frameon=False, handletextpad=0.1)


	ax[1].set_xlabel('Period [days]')
	# if per_resid.name == 'Bayesian Lomb-Scargle':
	# 	ax[1].set_ylabel('Normalized  P($f$ |D,I)')
	# 	m = max(per_resid.power)
	# 	ax[1].set_ylim([0,1.05])
	# else:
	ax[1].set_ylabel('Power')

	ax[1].semilogx(1./per_resid.freq, per_resid.power, 'k')
	ax[1].axhline(per_resid.powerFAP_1, c='k', lw=1, ls='--', dashes=(4,3), alpha=0.6, label=r'1\% FAP')
	ax[1].set_xlim([0.5, 20])
	ax[1].set_ylim([0, 20])
	ax[1].legend(frameon=False, handletextpad=0.1)
	# ax[1].xaxis.set_major_formatter(ScalarFormatter())

	ax[1].set_xticklabels(['', '', r'$1$', r'$10$'])
	# for t in ax[1].get_xticklabels():
	# 	t.set_text(t.get_text().replace('.0', ''))
	# print [x.get_text() for x in ticklabels]
	# ax[1].set_xticklabels(ticklabels)
	# ax[1].set_xticklabels(['1'])
	# rstyle(ax, grid_toggle=True)

	
	fig.tight_layout()

	try:
		sf = sys.argv[1]
		fig.savefig(sf)
	except IndexError:
		pass
plt.show()

