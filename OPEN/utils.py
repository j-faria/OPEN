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
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse
import matplotlib.cm as cm


mjup2mearth  = 317.828



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
    """ Root mean square of array"""
    return sqrt(sum(x*x for x in array)/len(array))

def wrms(array, weights):
    """ Root mean square of array"""
    w = weights #/ sum(weights)
    # return sqrt(sum(x*x for x in array*w))
    return sqrt(sum(w*(array - np.average(array, weights=w))**2) / sum(w))
    
 
def stdout_write(msg):
	""" Print to stdout (without carriage return) and flush right away.
	Useful to print in the same line """
	stdout.write(msg)
	stdout.flush()



### Matplotlib advanced plot interaction stuff
def selectable_plot(system, **kwargs):
    from shell_colors import yellow, blue
    from .logger import clogger
    msg = blue('INFO: ') + 'Click on a point in the plot to remove it.'
    clogger.info(msg)
    msg = blue('    : ') + 'Press ENTER when you are finished'
    clogger.info(msg)

    indices_to_remove = []
    global times
    times = 0
    def onpick3(event):
        global times

        ind = event.ind
        i, x, y = ind[0], np.take(system.time, ind)[0], np.take(system.vrad, ind)[0]
        indices_to_remove.append(i)

        # if times == 0:
        #     print
        #     msg = blue('INFO: ')
        #     times += 1
        # else:
        msg = blue('    : ')

        msg += 'going to remove observation %d -> %8.2f, %8.2f' % (i+1, x, y)
        clogger.info(msg)
        # print 'onpick3 scatter:', ind, np.take(system.time, ind), np.take(system.vrad, ind)
        # ax.scatter(np.take(system.time, ind), np.take(system.vrad, ind), color='r')
        # fig.show()

    fig, ax = plt.subplots()
    col = ax.scatter(system.time, system.vrad, picker=True)
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('RV [%s]'%system.units)
    fig.canvas.mpl_connect('pick_event', onpick3)
    plt.show()
    raw_input('')

    return unique(indices_to_remove)

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



### get number of cpus
def get_number_cores():
    """ Returns the number of available cpus in the machine """
    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass
    # http://code.google.com/p/psutil/
    try:
        import psutil
        return psutil.NUM_CPUS
    except (ImportError, AttributeError):
        pass
    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass
    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass
    # if all else fails
    return None

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
    fig = plt.figure(figsize=(dim, dim))
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
            plt.plot(y, x, '.')

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


def triangle_plot_kde(xs, labels=None, extents=None, truths=None, truth_color="#4682b4",
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
    import scipy.stats

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
    fig = plt.figure(figsize=(dim, dim))
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                        wspace=whspace, hspace=whspace)

    if extents is None:
        extents = [[x.min(), x.max()] for x in xs]

    for i, x in enumerate(xs):
        # Plot the histograms.
        ax = fig.add_subplot(K, K, i * (K + 1) + 1)

        gkde = scipy.stats.gaussian_kde(x)
        ind = np.linspace(min(x), max(x), 100)
        
        n, b, p = ax.hist(x, bins=kwargs.get("bins", 50), range=extents[i],
                histtype="step", color=kwargs.get("color", "k"))
        bin_width = np.ediff1d(b)[0] 

        ax.plot(ind, gkde.evaluate(ind)*bin_width*len(x), 'r')

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
            plt.plot(y, x, '.')

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