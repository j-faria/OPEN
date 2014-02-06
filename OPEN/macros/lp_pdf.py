# -*- coding: utf-8 -*-
# 'run -i' this script:
magic = get_ipython().magic

import datetime
import glob
import os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from OPEN.periodograms import gls
from OPEN.utils import day2year, rms
from OPEN.tqdm import tqdm


all_files = glob.glob('/home/joao/phd/data/*mean*')
# print len(all_files)

with PdfPages('lp_data.pdf') as pdf:
	for filename in tqdm(all_files[0:10]):
		# print filename
		print 
		magic('read ' + filename + ' --skip=2 -d')

		N = len(default.time)
		if N < 5: 
			print "LESS THAN 5 MEASUREMENTS"
			continue # skip this star

		per = gls(default) # calculate periodogram
		per1 = gls(default, quantity='fwhm')
		per2 = gls(default, quantity='rhk')

		tspan = max(default.time) - min(default.time)
		rms_over_err = rms(default.vrad - default.vrad.mean()) / default.error.mean()

		plt.figure(figsize=(8,10))
		gs = gridspec.GridSpec(4, 2)

		# plot the data
		ax1 = plt.subplot(gs[0])
		ax1.ticklabel_format(useOffset=False)
		default.do_plot_obs(newFig=False, leg=False)
		plt.locator_params(nbins=4)
		
		# info
		ax2 = plt.subplot(gs[0,1])
		plt.axis('off')
		plt.text(0., 0.9, os.path.basename(filename))
		plt.text(0., 0.68, '# points  %d' % N)
		plt.text(0., 0.5, 'time span  %3.1f years' % day2year(tspan))
		plt.text(0., 0.3, 'RV rms / <err>   %3.2f' % rms_over_err)
		# plt.text(0.5, 0.5, '# meas.  %d' % N)

		# plot the periodogram
		ax3 = plt.subplot(gs[1, :])
		per._plot(doFAP=True, newFig=False, axes=ax3)

		# the other periogorams	
		try:
			ax4 = plt.subplot(gs[2, :])
			per1._plot(newFig=False, axes=ax4)
			plt.title('FWHM')

			ax5 = plt.subplot(gs[3, :])
			per2._plot(newFig=False, axes=ax5)
			plt.title('RHK')
		except ValueError:
			print 'SOME ERROR OCCURRED...'
			continue


		pdf.savefig()
		plt.close()		


	# We can also set the file's metadata via the PdfPages object:
	d = pdf.infodict()
	d['Title'] = 'LP-metal-poor-data'
	d['Author'] = u'João Faria'
	d['Subject'] = 'Subject'
	d['Keywords'] = 'LP metal poor data gls joao faria'
	d['CreationDate'] = datetime.datetime.today()
	d['ModDate'] = datetime.datetime.today()
