# 'run -i' this script:
magic = get_ipython().magic

import matplotlib.pyplot as plt
import numpy as np
from OPEN.periodograms import gls

plt.figure(figsize=(8,10))
# two planets
magic('read /home/joao/HD41248_Kepler-88-2_2planets.dat -d')
real_periods = [10.95416, 70.18]

ax1 = plt.subplot(311)
p1 = gls(default, ofac=20)
p1._plot(verts=real_periods, color='k', doFAP=True, newFig=False, axes=ax1)
plt.title('Both planets')
plt.xlabel('')

magic('read /home/joao/HD41248_Kepler-88-2_planetXL.dat -d')

ax1 = plt.subplot(312)
p1 = gls(default, ofac=20)
p1._plot(verts=[real_periods[1],], color='k', doFAP=True, newFig=False, axes=ax1)
plt.title('Bigger planet')
plt.xlabel('')


magic('read /home/joao/HD41248_Kepler-88-2_planetXS.dat -d')

ax1 = plt.subplot(313)
p1 = gls(default, ofac=20)
p1._plot(verts=[real_periods[0],], color='k', doFAP=True, newFig=False, axes=ax1)
plt.title('Smaller planet')
# plt.xlabel('')


plt.show()