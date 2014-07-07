# -*- coding: utf-8 -*-
"""
This example demonstrates the use of pyqtgraph's dock widget system.

The dockarea system allows the design of user interfaces which can be rearranged by
the user at runtime. Docks can be moved, resized, stacked, and torn out of the main
window. This is similar in principle to the docking system built into Qt, but 
offers a more deterministic dock placement API (in Qt it is very difficult to 
programatically generate complex dock arrangements). Additionally, Qt's docks are 
designed to be used as small panels around the outer edge of a window. Pyqtgraph's 
docks were created with the notion that the entire window (or any portion of it) 
would consist of dockable components.

"""


import sys
import os
from time import time
sys.path.append('/home/joao/Work/Python/pyqtgraph-0.9.8')
sys.path.append('/home/joao/Work/OPEN')
from OPEN.classes import BasicTimeSeries
from OPEN.periodograms import gls
from OPEN.ext.get_rvN import get_rvn
import pyqtgraph as pg
# import pyqtgraph.multiprocess as mp
from pyqtgraph.Qt import QtGui
# import pyqtgraph.console
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType

import numpy as np
import random
import tables as tb

from pyqtgraph.dockarea import *

def create_planet_widget(star, database, debug=True):

	this_dir = os.path.dirname(os.path.realpath(__file__))
	data_relative_path = os.path.join(this_dir, '..', '..', 'data')
	available_sampling_files = {'random': [None, 0],
	                            '14Her': [os.path.join(data_relative_path, '14Her.times'), 203],
								'51Peg': [os.path.join(data_relative_path, '51Peg.times'), 409],
	                           }
	sampling_options = {'1': 'random', '2': "14Her", '3': "51Peg"}


	class PlanetCreator:
		def __init__(self, sampling_file, nobs):
			self.sampling_file = sampling_file
			self.nobs = nobs

			self.nplanets = 0
			self.periods = None
			self.eccentricities = None
			self.semi_amplitudes = None
			self.omegas = None
			self.times_periastron = None
			self.types_noise = []

			self.update_sampling(sampling_file, nobs)
			self.update_noise([], force=True)
			self.get_rv()

		def get_rv(self):
			self.vel_full = np.zeros_like(self.times_full)  # rv at oversampled times
			self.vel_full_each = np.zeros_like(self.times_full)  # rv at oversampled times for each planet
			self.vel_all = np.zeros((len(self.times), self.nplanets))  # rv at sample times for each planet
			self.vel_each = np.zeros_like(self.times)  # rv at sample times for each planet (temporary)
			self.vel_total = np.zeros_like(self.times)  # rv at sample times

			# if self.nplanets == 0:
			# 	self.get_noise(self.types_noise)
				# return

			for planet in range(self.nplanets):
				if self.periods is None: 
					P = np.random.rand()*998. + 2.  # random period, U(2, 998)
				else: 
					P = self.periods[planet]  # user-provided period

				if self.eccentricities is None:
					ecc = np.random.beta(0.867, 3.03)  # from Kipping 2013
				else:
					ecc = self.eccentricities[planet]  # user-provided eccentricity

				if self.semi_amplitudes is None:
					K = np.random.rand()*50. + 2.
				else:
					K = self.semi_amplitudes[planet]  # user-provided semi-amplitude

				if self.omegas is None:
					omega = np.random.rand()*2.*np.pi 
				else:
					omega = self.omegas[planet]  # user-provided omega

				if self.times_periastron is None:
					# random starting at around 2010-01-01, ending 1000 days later
					time_periastron = np.random.randint(2455197, 2456197) 
				else:
					time_periastron = self.times_periastron[planet]  # user-provided time of periastron

				# log
				# temp.write(output % (planet+1, P, K, ecc, omega, 2452000))

				# get RVs for this planet
				get_rvn(self.times, P, K, ecc, omega, time_periastron, 0., self.vel_each)

				# get RVs for this planet at oversampled times
				get_rvn(self.times_full, P, K, ecc, omega, time_periastron, 0., self.vel_full_each)

				# store
				self.vel_all[:,planet] = self.vel_each
				self.vel_total += self.vel_each
				self.vel_full += self.vel_full_each

			# self.get_noise(self.types_noise)
			self.vel_total += self.noise  # add noise


		def get_noise(self, types_noise):
			"""
			Return the noise vector to be added to the RVs
			"""
			noise_vector = np.zeros_like(self.times)

			for noise in types_noise:
				if noise['type'] == 'white':
					noise_vector += np.random.randn(len(self.times)) * noise['std']

			self.noise = noise_vector

		def ouput_to_file(self, filename):
			print filename

		def update_planets(self, nplanets):
			self.nplanets = nplanets
			self.get_rv()

		def update_noise(self, noise_dict, force=False):
			if (noise_dict == self.types_noise) and not force:
				return
			else:
				self.types_noise = noise_dict
				self.get_noise(noise_dict)

		def update_sampling(self, sampling_file, nobs):
			self.sampling_file = sampling_file
			self.nobs = nobs

			if sampling_file is not None:  # the sampling is set in the file
				times_sampled_from_file = np.loadtxt(sampling_file, usecols=(0,))
				print len(times_sampled_from_file)
				self.times = times_sampled_from_file
				self.times_full = np.linspace(min(times_sampled_from_file), max(times_sampled_from_file), 1000)
			else:  # the sampling is random
				self.times_full = np.linspace(2449460, 2452860, 1000)
				# sample (N points) randomly from the data to produce unevenly sampled time series
				rand_args = [i for i in sorted(random.sample(xrange(len(self.times_full)), nobs))]
				self.times = self.times_full[rand_args]



	PC = PlanetCreator(None, 10)
	TimeSer = BasicTimeSeries()
	TimeSer.time = PC.times
	TimeSer.vrad = PC.vel_total
	TimeSer.error = PC.noise
	per = gls(TimeSer)



	app = QtGui.QApplication.instance()  # checks if QApplication already exists
	if not app:  # create QApplication if it doesnt exist 
		app = QtGui.QApplication([])
	win = QtGui.QMainWindow()
	area = DockArea()
	win.setCentralWidget(area)
	win.resize(1250,550)
	win.setWindowTitle('OPEN')

	## Create docks, place them into the window one at a time.
	## Note that size arguments are only a suggestion; docks will still have to
	## fill the entire dock area and obey the limits of their internal widgets.
	d1 = Dock("Info (d1)", size=(300, 100))     ## give this dock the minimum possible size
	d3 = Dock("Options (d3)", size=(300,400))
	d4 = Dock("per (d4)", size=(500,100))
	d6 = Dock("RV (d6)", size=(500,300))

	area.addDock(d1, 'left')      
	area.addDock(d3, 'bottom', d1)## place d3 at bottom edge of d1
	area.addDock(d4, 'bottom')     ## place d4 at right edge of dock area
	area.addDock(d6, 'right')   ## place d5 at top edge of d4

	## Test ability to move docks programatically after they have been placed
	area.moveDock(d4, 'bottom', d6)   ## move d6 to stack on top of d4


	## this group includes a menu allowing the user to add new parameters into its child list
	class PlanetParametersGroup(pTypes.GroupParameter):
		def __init__(self, **opts):
			opts['type'] = 'group'
			opts['addText'] = "Add"
			opts['addList'] = ['planet']
			pTypes.GroupParameter.__init__(self, **opts)

		def addNew(self, typ):
			planet_group = [
            				{'name': 'period [d]', 'type': 'float', 'value': 10},
            				{'name': 'K [m/s]', 'type': 'float', 'value': 10},
            				{'name': 'eccentricity', 'type': 'float', 'value': 0.0, 'step': 0.01},
            				{'name': 'w [rad]', 'type': 'float', 'value': 0.0, 'step': 0.01},
            				{'name': 'Tp [d]', 'type': 'float', 'value': 2455197, 'step': 10},
            				{'name': 'Randomize', 'type': 'action'},
            			   ]
			self.addChild(dict(name="planet %d" % (len(self.childs)+1), 
				               type='group', children=planet_group, removable=True, renamable=False))

	## this group includes a menu allowing the user to add new parameters into its child list
	class NoiseParametersGroup(pTypes.GroupParameter):
		def __init__(self, **opts):
			opts['type'] = 'group'
			opts['addText'] = "Add"
			opts['addList'] = ['white', 'correlated']
			pTypes.GroupParameter.__init__(self, **opts)

		def addNew(self, typ):
			print typ
			white_noise_group = [{'name': 'variance', 'type': 'float', 'value': 1, 'step': 0.1}]

			if typ == 'white':
				self.addChild(dict(name="white noise %d" % (len(self.childs)+1), 
					               type='group', children=white_noise_group, removable=True))
			elif typ == 'correlated':
				pass


	## User-changeable parameters
	params = [
        	# number of points
        	{'name': 'Number of points', 'type': 'int', 'value': PC.nobs},
        	# parameters for planets
        	PlanetParametersGroup(name="Planet parameters", children=[]),
        	# choose sampling times
        	{'name': 'Sampling times', 'type': 'group', 'children': [
        		{'name': 'Sampling', 'type': 'list', 'values': {"Random": 1, "14Her": 2, "51Peg": 3}, 'value': 1},
        		{'name': 'Start (dw)', 'type': 'float', 'value': 2449460},
        		{'name': 'End (dw)', 'type': 'float', 'value': 2452860},
        		]},
        	# parameters for noise
        	NoiseParametersGroup(name='Noise Parameters', children=[]),
        	# output to file
        	{'name': 'Save', 'type': 'group', 'children': [
        		{'name': 'Filename', 'type': 'str', 'value': "example.rdb"},
        		{'name': 'Save', 'type': 'action'},
        		]},
        	]


	## Add widgets into each dock

	## first dock gets info
	w1 = pg.LayoutWidget()
	label_text = \
	""" OPEN -- %s	B: %.2f	  V: %.2f 
	%s	FeH		logg
	%dK	%.2f+-%.2f	%.2f+-%.2f
	Data: %d measurements	%s
	""" 

	teff = feh = feh_error = logg = logg_error = 0

	label = QtGui.QLabel(label_text % (star, 0, 0,
		                               'G2', 
		                               teff, feh, feh_error, logg, logg_error,
		                               100, 'HARPS')
	                    )


	## Create tree of Parameter objects
	p = Parameter.create(name='params', type='group', children=params)

	def clear_plots():
		w4.clear()
		w6.clear()

	def do_replot():

		### RVs ============================
		time = PC.times
		rv = PC.vel_total
		errors = PC.noise

		### calculate periodogram
		TimeSer.time = PC.times
		TimeSer.vrad = PC.vel_total
		TimeSer.error = PC.noise
		per = gls(TimeSer)
		print per.power

		# teff = feh = feh_error = logg = logg_error = 0
		# l = w1.children()[1]
		# l.setText(label_text % (star, 0, 0,
		# 	                     'G2', 
		# 	                     teff, feh, feh_error, logg, logg_error,
		# 	                     100, 'HARPS')
		# 		   )

		# w1.addWidget(label, row=0, col=0)

		## dock 4, periodogram
		w4.plot(1./per.freq, per.power)
		w4.autoRange()


		## dock 6, data
		w6.plot(PC.times_full, PC.vel_full)
		w6.plot(time, rv, pen=None, symbolBrush=(255, 0, 0), symbolSize=10)
		errItem = pg.ErrorBarItem(x=time, y=rv, height=errors, beam=0.5)
		w6.addItem(errItem)
		w6.autoRange()


	## If anything changes in the tree
	def change(param, changes):
		nobs = param['Number of points']
		signals = param.children()[1]
		noises = param.children()[3]

		nplanets = len(signals.children())
		# nnoises = len(noises.children())
		# print signals.children()
		print changes

		noise_dict = []
		for n in noises:
			if 'white noise' in n.name():
				noise_dict.append({'type': 'white', 'std': np.sqrt(n['variance'])})
			else:
				pass
		PC.update_noise(noise_dict)

		# status of "Sampling times" group
		sampling = str(param.children()[2].children()[0].value())
		sampling_file = available_sampling_files[sampling_options[sampling]][0]
		sampling_file_npoints = available_sampling_files[sampling_options[sampling]][1]
		if sampling_file_npoints > 0:
			param['Number of points'] = nobs = sampling_file_npoints

		# planet signals
		if nplanets > 0:
			PC.periods = []
			PC.eccentricities = []
			PC.semi_amplitudes = []
			PC.omegas = []
			PC.times_periastron = []

		for s in signals:
			# print s.children()
			PC.periods.append(s['period [d]'])
			PC.eccentricities.append(s['eccentricity'])
			PC.semi_amplitudes.append(s['K [m/s]'])
			PC.omegas.append(s['w [rad]'])
			PC.times_periastron.append(s['Tp [d]'])
			# print s.children()
		PC.update_planets(nplanets)


		save = param.children()[4].children()[0]
		filename = save.value()


		# print("tree changes:")
		for param, change, data in changes:
			path = p.childPath(param)
			if path is not None:
				childName = '.'.join(path)
			else:
				childName = param.name()
			# print('  parameter: %s' % childName)
			# print('  change:    %s' % change)
			# print('  data:      %s' % str(data))
			# print('  ----------')

			if childName == 'Save.Save' and change == 'activated':
				PC.ouput_to_file(filename)


		if (nobs != PC.nobs) or (sampling_file != PC.sampling_file):
			PC.update_sampling(sampling_file, nobs)
			PC.update_noise(noise_dict, force=True)
			PC.update_planets(nplanets)

		clear_plots()
		do_replot()
    
	p.sigTreeStateChanged.connect(change)



	## Create ParameterTree widget
	tree = ParameterTree()
	tree.setParameters(p, showTop=False)


	w1.addWidget(label, row=0, col=0)
	# w1.addWidget(saveBtn, row=1, col=0)
	# w1.addWidget(restoreBtn, row=1, col=1)
	# w1.addWidget(tree, row=0, col=2, rowspan=1)
	d1.addWidget(w1)

	d3.addWidget(tree)
	# state = None
	# def save():
	#     global state
	#     state = area.saveState()
	#     restoreBtn.setEnabled(True)
	#     print state
	# def load():
	#     global state
	#     area.restoreState(state)
	# saveBtn.clicked.connect(save)
	# restoreBtn.clicked.connect(load)


	## Hide title bar on dock 3, the radial velocity data
	# d3.hideTitleBar()
	# w3 = pg.PlotWidget(title="RV data", name="PlotRV")
	# w3.plot(x, y, pen=None, symbol='o', symbolSize=5)
	# errItem = pg.ErrorBarItem(x=x, y=y, height=err, beam=0.5)
	# w3.addItem(errItem)
	# d3.addWidget(w3)


	w4 = pg.PlotWidget(title="periodogram ")
	# w4.plot(np.random.normal(size=100))
	w4.plot(1./per.freq, per.power)
	w4.plotItem.setLogMode(True, False)
	d4.addWidget(w4)

	### RVs ============================
	time = PC.times
	rv = PC.vel_total
	errors = PC.noise

	w6 = pg.PlotWidget(name="perRV")
	w6.plot(PC.times_full, PC.vel_full)
	w6.plot(time, rv, pen=None, symbolBrush=(255,0,0), symbolSize=10)
	errItem = pg.ErrorBarItem(x=time, y=rv, height=errors, beam=0.5)
	w6.addItem(errItem)

	# w6.plotItem.setLogMode(True, False)
	#cross hair
	# vLine = pg.InfiniteLine(angle=90, movable=False)
	# hLine = pg.InfiniteLine(angle=0, movable=False)
	# w6.addItem(vLine, ignoreBounds=True)
	# w6.addItem(hLine, ignoreBounds=True)
	# vb = w6.plotItem.vb

	# def mouseMoved(evt):
	# 	pos = evt[0]  ## using signal proxy turns original arguments into a tuple
	# 	if w6.sceneBoundingRect().contains(pos):
	# 		mousePoint = vb.mapSceneToView(pos)
	# 		vLine.setPos(mousePoint.x())
	# 		hLine.setPos(mousePoint.y())

	# proxy = pg.SignalProxy(w6.scene().sigMouseMoved, rateLimit=60, slot=mouseMoved)
	#p1.scene().sigMouseMoved.connect(mouseMoved)
	d6.addWidget(w6)


	### periodogram RHK
	# w8 = pg.PlotWidget()
	# w8.plot(1./freq, power_rhk, title="Periodogram FWHM")
	# w8.plotItem.setLogMode(True, False)
	# w8.setXLink('perRV')  ## test linking by name
	# d8.addWidget(w8)

	win.show()
	QtGui.QApplication.instance().exec_()

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    create_planet_widget('HD41248', None)

    # if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        # QtGui.QApplication.instance().exec_()

## Running from inside OPEN
if __name__ == 'OPEN.plugins.create_GUI.run':
	create_planet_widget('HD190984', None)
