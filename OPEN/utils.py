# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
"""
Utility functions or snippets for all sorts of things
"""

from math import sqrt
from sys import stdout


## this code from 
##   http://code.activestate.com/recipes/502263-yet-another-unique-function/
## atributed to Paul Rubin
################################################################################
def unique(seq, keepstr=True):
    """ Return unique elements from enumerables (list, tuple, str)"""
    t = type(seq)
    if t==str:
        t = (list, ''.join)[bool(keepstr)]
    seen = []
    return t(c for c in seq if not (c in seen or seen.append(c)))



def day2year(day):
    """ Convert days to years """
    return day*0.00273791
    
    
def rms(array):
    return sqrt(sum(n*n for n in array)/len(array))
    
    
def stdout_write(msg):
	""" Print to stdout (without carriage return) and flush right away.
	Useful to print in the same line """
	stdout.write(msg)
	stdout.flush()



### Matplotlib advanced plot interaction stuff
from matplotlib.widgets import RectangleSelector
import numpy as np
import matplotlib.pyplot as plt

def line_select_callback(eclick, erelease):
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print ("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    print (" The button you used were: %s %s" % (eclick.button, erelease.button))

def toggle_selector(event):
    print (' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print (' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print (' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)

def selectable_plot(*args, **kwargs):
    fig, current_ax = plt.subplots()                    # make a new plotingrange

    plt.plot(*args, **kwargs)  # plot something

    # drawtype is 'box' or 'line' or 'none'
    toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1,3], # don't use middle button
                                           minspanx=5, minspany=5,
                                           spancoords='pixels')
    plt.connect('key_press_event', toggle_selector)
    plt.show()