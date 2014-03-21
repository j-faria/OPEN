# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
"""
Utility functions or snippets for all sorts of things
"""

from math import sqrt, ceil
from sys import stdout
import os
import subprocess

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse
import matplotlib.cm as cm

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

def julian_day_to_date(J):
    """ Returns the date corresponding to a julian day number"""
    J = int(ceil(J))  # i think we need to do this...
    y=4716; v=3; j=1401; u=5; m=2; s=153; n=12; w=2; r=4; B=274277; p=1461; C=-38
    if len(str(int(J))) < 7: J = J+2400000
    f = J + j + (((4 * J + B)/146097) * 3)/4 + C
    e = r * f + v
    g = (e % p)/r
    h = u * g + w
    D = (h % s)/u + 1
    M = ((h/s + m) % n) + 1
    Y = e/p - y + (n + m - M)/n
    return (Y, M, D)


### Yorbit related I/O
def write_yorbit_macro(system):
    with open('./OPEN/yorbit.macro.template') as f:
        template = f.read()

    filenames = []
    for i,f in enumerate(system.provenance.keys()):
        cmd = 'cp ' + f + ' ' + '/home/joao/yorbit/data/rv_files/rv_file%d.rdb' % (i+1)
        filenames.append('rv_file%d.rdb' % (i+1))
        os.system(cmd)
        #print cmd
    #print filenames

    model = ''
    try:
        for key, val in system.model.iteritems():
            if val >0: model += '%s%d ' % (key, val)
    except AttributeError:
        model = raw_input('Which model? ')

    with open('/home/joao/yorbit/work/macros/macro1.macro', 'w') as f:
        print >>f, template % (' '.join(filenames), '', model)

    cmd = ['/home/joao/Software/yorbit/yorbit-1.4.6/src/batch/batch_macro.i']
    cmd += ['/home/joao/Software/yorbit/yorbit-1.4.6/src']
    cmd += ['macro1.macro']
    # print cmd
    with open('yorbit.log', "w") as outfile:
        subprocess.call(cmd, stdout=outfile)




def ask_yes_no(msg, default):
    # raw_input returns the empty string for "enter"
    yes = set(['yes','y','ye','Y'])
    no = set(['no','n','NO','N'])
    print msg,
    choice = raw_input().lower()
    if not choice:
        return default
    if choice in yes:
        return True
    elif choice in no:
        return False
    else:
        stdout.write("Please respond with 'yes' or 'no'")
        return False


def notnan(v):
    return v[~np.isnan(v)]


def triangle_plot(xs, labels=None, extents=None, truths=None, truth_color="#4682b4",
                  scale_hist=False, quantiles=[], **kwargs):
    """
    Make a *sick* corner plot showing the projections of a set of samples
    drawn in a multi-dimensional space.

    :param xs: ``(nsamples, ndim)``
        The samples. This should be a 1- or 2-dimensional array. For a 1-D
        array this results in a simple histogram. For a 2-D array, the zeroth
        axis is the list of samples and the next axis are the dimensions of
        the space.

    :param labels: ``ndim`` (optional)
        A list of names for the dimensions.

    :param truths: ``ndim`` (optional)
        A list of reference values to indicate on the plots.

    :param truth_color: (optional)
        A ``matplotlib`` style color for the ``truths`` makers.

    :param quantiles: (optional)
        A list of fractional quantiles to show on the 1-D histograms as
        vertical dashed lines.

    :param scale_hist: (optional)
        Should the 1-D histograms be scaled in such a way that the zero line
        is visible?

    """

    # Deal with 1D sample lists.
    xs = np.atleast_1d(xs)
    if len(xs.shape) == 1:
        xs = np.atleast_2d(xs)
    else:
        assert len(xs.shape) == 2, "The input sample array must be 1- or 2-D."
        xs = xs.T
    assert xs.shape[0] <= xs.shape[1], "I don't believe that you want more " \
                                       "dimensions than samples!"

    K = len(xs)
    factor = 2.0           # size of one side of one panel
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.05 * factor  # size of top/right margin
    whspace = 0.05         # w/hspace size
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim = lbdim + plotdim + trdim
    fig = pl.figure(figsize=(dim, dim))
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                        wspace=whspace, hspace=whspace)

    if extents is None:
        extents = [[x.min(), x.max()] for x in xs]

    for i, x in enumerate(xs):
        # Plot the histograms.
        ax = fig.add_subplot(K, K, i * (K + 1) + 1)
        n, b, p = ax.hist(x, bins=kwargs.get("bins", 50), range=extents[i],
                histtype="step", color=kwargs.get("color", "k"))
        if truths is not None:
            ax.axvline(truths[i], color=truth_color)

        # Plot quantiles if wanted.
        if len(quantiles) > 0:
            xsorted = sorted(x)
            for q in quantiles:
                ax.axvline(xsorted[int(q * len(xsorted))], ls="dashed",
                           color=kwargs.get("color", "k"))

        # Set up the axes.
        ax.set_xlim(extents[i])
        if scale_hist:
            maxn = np.max(n)
            ax.set_ylim(-0.1 * maxn, 1.1 * maxn)
        else:
            ax.set_ylim(0, 1.1 * np.max(n))
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(MaxNLocator(5))

        # Not so DRY.
        if i < K - 1:
            ax.set_xticklabels([])
        else:
            [l.set_rotation(45) for l in ax.get_xticklabels()]
            if labels is not None:
                ax.set_xlabel(labels[i])
                ax.xaxis.set_label_coords(0.5, -0.3)

        for j, y in enumerate(xs[:i]):
            ax = fig.add_subplot(K, K, (i * K + j) + 1)
            pl.plot(y, x, '.')

            if truths is not None:
                ax.plot(truths[j], truths[i], "s", color=truth_color)
                ax.axvline(truths[j], color=truth_color)
                ax.axhline(truths[i], color=truth_color)

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_xticklabels()]
                if labels is not None:
                    ax.set_xlabel(labels[j])
                    ax.xaxis.set_label_coords(0.5, -0.3)

            if j > 0:
                ax.set_yticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_yticklabels()]
                if labels is not None:
                    ax.set_ylabel(labels[i])
                    ax.yaxis.set_label_coords(-0.3, 0.5)

    return fig