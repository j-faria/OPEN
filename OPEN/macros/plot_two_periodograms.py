# In order to plot the periodogram of the RV and the FWHM together 'run -i' this script:
magic = get_ipython().magic

import matplotlib.pyplot as plt

# magic('per obs --bayes')  # to calculate GLS for RV
# magic('per fwhm --bayes')  # to calculate GLS for RV

magic('per obs')  # to calculate GLS for RV
magic('per fwhm')  # to calculate GLS for RV

per_rv = default.per
per_fwhm = default.fwhm_per


plt.figure()
gs = gridspec.GridSpec(2, 1)
ax = [plt.subplot(gs[0]), plt.subplot(gs[1])]
# plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.15, hspace=0.4)

if per_rv.name == 'Bayesian Lomb-Scargle':
	ax[0].set_ylabel(r'Normalized  P($f$ |D,I)')
	m = max(per_rv.power)
	ax[0].set_ylim([0,1.05])
else:
	ax[0].set_ylabel('Normalized Power')
	m = 1.

ax[0].semilogx(1./per_rv.freq, per_rv.power / m, 'k', label='RV')
ax[0].legend()

ax[1].set_xlabel('Period [days]')
if per_fwhm.name == 'Bayesian Lomb-Scargle':
	ax[1].set_ylabel('Normalized  P($f$ |D,I)')
	m = max(per_fwhm.power)
	ax[1].set_ylim([0,1.05])
else:
	ax[1].set_ylabel('Normalized Power')
	m = 1.

ax[1].semilogx(1./per_fwhm.freq, per_fwhm.power / m, 'k', label='FWHM')
ax[1].legend()

rstyle(ax, grid_toggle=True)

plt.show()