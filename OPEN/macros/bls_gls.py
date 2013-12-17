# In order to plot the GLS and BLS together 'run -i' this script:
magic = get_ipython().magic

magic('per obs')  # to calculate GLS

from OPEN.periodograms import bls
default.per2 = bls(default)  # calculate BLS and store it in the system

# for normalization
a1 = default.per2.power.max()
a2 = default.per.power.max()

from matplotlib.pylab import semilogx
semilogx(1./default.per.freq, default.per.power, 'b-', label='gls')
semilogx(1./default.per2.freq, default.per2.power/a1*a2, 'r-', label='gls')

