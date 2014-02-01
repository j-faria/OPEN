# 'run -i' this script:
magic = get_ipython().magic
from OPEN.periodograms import gls, bls
import matplotlib.pyplot as plt

magic('read /home/joao/yorbit/data/rv_files/HD41248_omc_d1_harps_mean_set1.rdb --skip=2 -v -d')
per1a = gls(default, ofac=12)
per1b = bls(default)

magic('read /home/joao/yorbit/data/rv_files/HD41248_omc_d1_harps_mean_set2.rdb /home/joao/yorbit/data/rv_files/HD41248_omc_d1_harps_mean_set3.rdb --skip=2 -v -d')
per2 = gls(default, ofac=12)
per2a = bls(default)

# magic('read /home/joao/yorbit/data/rv_files/HD41248_omc_d1_harps_mean_set3.rdb --skip=2 -v -d')
# per3 = gls(default, ofac=12)


# for normalization
a1 = per1b.power.max()
a2 = per1a.power.max()

plt.figure()
plt.semilogx(1./per1a.freq, per1a.power, 'b-', label='set1-gls')
plt.semilogx(1./per1b.freq, per1b.power, 'b:', label='set1-bls')
plt.semilogx(1./per2a.freq, per2a.power, 'r-', label='set2-gls')
plt.semilogx(1./per2b.freq, per2b.power/a1*a2, 'r:', label='set2-bls')
# plt.semilogx(1./per3.freq, per3.power, 'g-', label='set3')
plt.legend()
plt.show()