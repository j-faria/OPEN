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
from time import time
sys.path.append('/home/joao/Work/Python/pyqtgraph-0.9.8')
sys.path.append('/home/joao/Work/OPEN')
from OPEN.classes import rvSeries
from OPEN.periodograms import gls
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.console
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType

import numpy as np
import tables as tb

from pyqtgraph.dockarea import *

def metalpoor_widget(star, database, debug=False):
	# data from database
	database_file = '/home/joao/phd/data/metal-poor.h5'
	h5file = tb.open_file(database_file, mode="r")
	available_stars = [n._v_name for n in h5file.iter_nodes('/')]

	def get_data(h5file, star):
		t0 = time()

		d = getattr(h5file.root, star)
		attrs = d._v_attrs
		# global x, y, err 
		x = np.array([row['time'] for row in d.data])
		y = np.array([row['vrad'] for row in d.data])
		err = np.array([row['svrad'] for row in d.data])

		freq = np.array([row['freq'] for row in d.analysis])
		power_rv = np.array([row['power_RV'] for row in d.analysis])
		power_fwhm = np.array([row['power_FWHM'] for row in d.analysis])
		power_bis = np.array([row['power_BIS'] for row in d.analysis])
		power_rhk = np.array([row['power_RHK'] for row in d.analysis])

		try:
			bfreq = np.array([row for row in d.bfreq])
			bpower_rv = np.array([row for row in d.bpower_rv])
		except tb.NoSuchNodeError:
			bfreq = None
			bpower_rv = None

		if debug: print 'array creation took: %.3f s' % (time()-t0)

		return d, attrs, \
		       x,y,err, \
		       freq, power_rv, power_fwhm, power_bis, power_rhk,\
		       bfreq, bpower_rv

	d, attrs, \
		x,y,err, \
		freq,power_rv,power_fwhm,power_bis,power_rhk, \
		bfreq, bpower_rv = get_data(h5file, star)

	app=QtGui.QApplication.instance() # checks if QApplication already exists
	if not app: # create QApplication if it doesnt exist 
		app = QtGui.QApplication([])
	win = QtGui.QMainWindow()
	area = DockArea()
	win.setCentralWidget(area)
	win.resize(1250,550)
	win.setWindowTitle('OPEN')

	## Create docks, place them into the window one at a time.
	## Note that size arguments are only a suggestion; docks will still have to
	## fill the entire dock area and obey the limits of their internal widgets.
	d1 = Dock("Info", size=(1, 1))     ## give this dock the minimum possible size
	d2 = Dock("Console", size=(500,100))
	d3 = Dock("Dock3", size=(500,400))
	d4 = Dock("Info RV", size=(500,300))
	d5 = Dock("Periodogram - BIS", size=(500,200))
	d6 = Dock("Periodogram - RV", size=(500,300))
	d7 = Dock("Periodogram - FWHM", size=(500,200))
	d8 = Dock("Periodogram - RHK", size=(500,200))
	area.addDock(d1, 'left')      ## place d1 at left edge of dock area (it will fill the whole space since there are no other docks yet)
	area.addDock(d2, 'right')     ## place d2 at right edge of dock area
	area.addDock(d3, 'bottom', d1)## place d3 at bottom edge of d1
	area.addDock(d4, 'right')     ## place d4 at right edge of dock area
	area.addDock(d5, 'left', d1)  ## place d5 at left edge of d1
	area.addDock(d6, 'top', d4)   ## place d5 at top edge of d4
	area.addDock(d7, 'top', d2)
	area.addDock(d8, 'top', d2)

	## Test ability to move docks programatically after they have been placed
	area.moveDock(d4, 'top', d2)     ## move d4 to top edge of d2
	area.moveDock(d6, 'above', d4)   ## move d6 to stack on top of d4
	area.moveDock(d5, 'top', d2)     ## move d5 to top edge of d2
	area.moveDock(d7, 'above', d5)   ## move 
	area.moveDock(d8, 'above', d5)   ## move 


	## User-changeable parameters
	params = [
        	{'name': 'Star:', 'type': 'list', 'values': available_stars, 'value': star},
        	# {'name': 'Load', 'type': 'action'},
		]

	## Add widgets into each dock

	## first dock gets save/restore buttons
	w1 = pg.LayoutWidget()
	label_text = \
	""" OPEN -- metalpoor sample 

Star: %s	B: %.2f	  V: %.2f
	
	%s	FeH		logg
	%dK	%.2f+-%.2f	%.2f+-%.2f

Data: %d measurements	%s
""" 

	try: 
		teff, feh, feh_error, logg, logg_error = attrs.teff, attrs.feh, attrs.feh_error, attrs.logg, attrs.logg_error
	except AttributeError:
		teff = feh = feh_error = logg = logg_error = 0

	label = QtGui.QLabel(label_text % (star, attrs.bmag, attrs.vmag,
		                               attrs.spect_type, 
		                               teff, feh, feh_error, logg, logg_error,
		                               len(x), 'HARPS')
	                    )


	saveBtn = QtGui.QPushButton('Save dock')
	restoreBtn = QtGui.QPushButton('Restore dock')
	restoreBtn.setEnabled(False)

	## Create tree of Parameter objects
	p = Parameter.create(name='params', type='group', children=params)

	def clear_plots():
		w3.clear()
		w5.clear()
		w6.clear()
		w7.clear()
		w8.clear()

	def do_replot(star):
		d, attrs, \
			x,y,err, \
			freq,power_rv,power_fwhm,power_bis,power_rhk, \
			bfreq, bpower_rv = get_data(h5file, star)

		# w1.removeItem(label)
		try: 
			teff, feh, feh_error, logg, logg_error = attrs.teff, attrs.feh, attrs.feh_error, attrs.logg, attrs.logg_error
		except AttributeError:
			teff = feh = feh_error = logg = logg_error = 0
		l = w1.children()[1]
		l.setText(label_text % (star, attrs.bmag, attrs.vmag,
			                     attrs.spect_type, 
			                     teff, feh, feh_error, logg, logg_error,
			                     len(x), 'HARPS')
				   )

		# w1.addWidget(label, row=0, col=0)

		## dock 3, the radial velocity data
		w3.plot(x, y, pen=None, symbol='o', symbolSize=5)
		errItem = pg.ErrorBarItem(x=x, y=y, height=err, beam=0.5)
		w3.addItem(errItem)
		w3.autoRange()

		## dock 5, periodogram BIS
		w5.plot(1./freq, power_bis, title="Periodogram BIS")

		## dock 6, periodogram RVs
		w6.plot(1./freq, power_rv, title="Periodogram")
		if bfreq is not None:
			w6.plot(1./bfreq, bpower_rv*max(power_rv), pen='r')
		w6.autoRange()
		w6.addItem(vLine, ignoreBounds=True)
		w6.addItem(hLine, ignoreBounds=True)
		
		## dock 7, periodogram FWHM
		w7.plot(1./freq, power_fwhm, title="Periodogram FWHM")

		## dock 8, periodogram RHK
		w8.plot(1./freq, power_rhk, title="Periodogram FWHM")


	## If anything changes in the tree
	def change(param, changes):
		param, change, data = changes[0]
		if debug: 
			print 'Searching for star',

			for param, change, data in changes:
				print data,
			print '...',

		clear_plots()
		do_replot(data)
    
	p.sigTreeStateChanged.connect(change)



	## Create ParameterTree widget
	tree = ParameterTree()
	tree.setParameters(p, showTop=False)


	w1.addWidget(label, row=0, col=0)
	# w1.addWidget(saveBtn, row=1, col=0)
	# w1.addWidget(restoreBtn, row=1, col=1)
	w1.addWidget(tree, row=0, col=2, rowspan=1)
	d1.addWidget(w1)
	state = None
	def save():
	    global state
	    state = area.saveState()
	    restoreBtn.setEnabled(True)
	    print state
	def load():
	    global state
	    area.restoreState(state)
	saveBtn.clicked.connect(save)
	restoreBtn.clicked.connect(load)

	## Second dock gets the console widget
	w2 = pg.console.ConsoleWidget()
	d2.addWidget(w2)

	## Hide title bar on dock 3, the radial velocity data
	d3.hideTitleBar()
	w3 = pg.PlotWidget(title="RV data", name="PlotRV")
	w3.plot(x, y, pen=None, symbol='o', symbolSize=5)
	errItem = pg.ErrorBarItem(x=x, y=y, height=err, beam=0.5)
	w3.addItem(errItem)
	d3.addWidget(w3)


	# w4 = pg.PlotWidget(title="periodogram ")
	# w4.plot(np.random.normal(size=100))
	# d4.addWidget(w4)

	### periodogram RVs ============================
	w6 = pg.PlotWidget(name="perRV")
	w6.plot(1./freq, power_rv, title="Periodogram")
	if bfreq is not None:
		w6.plot(1./bfreq, bpower_rv*max(power_rv), pen='r')
	w6.plotItem.setLogMode(True, False)
	#cross hair
	vLine = pg.InfiniteLine(angle=90, movable=False)
	hLine = pg.InfiniteLine(angle=0, movable=False)
	w6.addItem(vLine, ignoreBounds=True)
	w6.addItem(hLine, ignoreBounds=True)
	vb = w6.plotItem.vb

	def mouseMoved(evt):
		pos = evt[0]  ## using signal proxy turns original arguments into a tuple
		if w6.sceneBoundingRect().contains(pos):
			mousePoint = vb.mapSceneToView(pos)
			vLine.setPos(mousePoint.x())
			hLine.setPos(mousePoint.y())

	proxy = pg.SignalProxy(w6.scene().sigMouseMoved, rateLimit=60, slot=mouseMoved)
	#p1.scene().sigMouseMoved.connect(mouseMoved)
	d6.addWidget(w6)

	### periodogram BIS
	w5 = pg.PlotWidget()
	w5.plot(1./freq, power_bis, title="Periodogram BIS")
	w5.plotItem.setLogMode(True, False)
	w5.setXLink('perRV')  ## test linking by name
	d5.addWidget(w5)


	### periodogram FWHM
	w7 = pg.PlotWidget()
	w7.plot(1./freq, power_fwhm, title="Periodogram FWHM")
	w7.plotItem.setLogMode(True, False)
	w7.setXLink('perRV')  ## test linking by name
	d7.addWidget(w7)

	### periodogram RHK
	w8 = pg.PlotWidget()
	w8.plot(1./freq, power_rhk, title="Periodogram FWHM")
	w8.plotItem.setLogMode(True, False)
	w8.setXLink('perRV')  ## test linking by name
	d8.addWidget(w8)



	win.show()
	QtGui.QApplication.instance().exec_()
	h5file.close()
	# return d6


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    metalpoor_widget('HD41248', None)

    # if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        # QtGui.QApplication.instance().exec_()

## Running from inside OPEN
if __name__ == 'OPEN.plugins.metalpoor.run':
	metalpoor_widget('HD190984', None)
