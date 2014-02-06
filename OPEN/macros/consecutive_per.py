# 'run -i' this script:
magic = get_ipython().magic

import matplotlib.pyplot as plt
import numpy as np
from OPEN.periodograms import gls
from matplotlib.backends.backend_pdf import PdfPages

magic('read /home/joao/yorbit/data/rv_files/HD41248_harps_mean.rdb --skip=2 -v -d')

with PdfPages('consecutive_periodograms.pdf') as pdf:

	for i in reversed(range(155, len(default.time)+1)):
		per = gls(default)
		print 'gls %d points from beginning' % i
		per._plot(verts=[18, 25])
		plt.title('gls %d points from beginning' % i)
		plt.ylim([0, 30])
		plt.xlim([0.1, 1e4])
		pdf.savefig()
		plt.close()
		default.time = np.delete(default.time, -1)
		default.vrad = np.delete(default.vrad, -1)
		default.error = np.delete(default.error, -1)
