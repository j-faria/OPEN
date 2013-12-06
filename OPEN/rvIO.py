# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

# standard library imports
from fileinput import FileInput
from itertools import islice, chain
import os

# other imports
from numpy import loadtxt

# intra-package imports
from .logger import clogger, logging
from .utils import day2year, rms

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)


def read_rv(*filenames, **kwargs):
    """
    Read one or more files containing radial velocity measurements.
    
    Parameters
    ----------
    filenames: string
        One or more files to read.
    
    Returns
    -------
    t: array
        Times of observations.
    rv: array
        Radial velocity values.
    err: array
        Error in the radial velocity.
    dic: dict
        Dictionary with name and number of values from each file.
    """

    # set logging level
    clogger.setLevel(logging.VERBOSE) \
                if (kwargs.has_key('verbose') and kwargs['verbose']) \
                else clogger.setLevel(logging.INFO)

    # how many header lines to skip?
    if (kwargs.has_key('skip')): header_skip = int(kwargs['skip'])

    dic = {} # will hold how many values per file
    for filename in sorted(filenames):
        if os.path.isfile(filename) and os.access(filename, os.R_OK):
            # this file exists and is readable
            with rvfile(filename) as f:
                nlines = len(f.readuncommented())
            dic[filename] = [nlines, 0]
            clogger.info('Reading %d values from file %s' % (nlines, filename))
        else:
            # should raise an error or read from the other files?
            raise IOError("The file '%s' doesn't seem to exist" % filename)
            
    # black magic to build input from file list while skipping headers
    finput = [FileInput(f) for f in sorted(filenames)]
    iterables = [islice(f, header_skip, None) for f in finput]
    files = chain(*iterables)

    # read data
    try:
        t, rv, err = loadtxt(files, unpack=True, usecols=(0,1,2))
    except Exception as e:
        raise e
    
    # verbose stats about data
    stats = None
    if (kwargs.has_key('verbose') and kwargs['verbose']):
        tspan = max(t) - min(t)
        stats = '\n'
        stats += "Total timespan : %f days = %f years\n" % (tspan, day2year(tspan))
        stats += "RV rms [m/s] : %f\n etc etc" % rms(rv)
    
    clogger.verbose(stats)
    
    return t, rv, err, dic
    
    
    
class rvfile(file):
    """
    Subclass of Python's File class that implements specific methods
    for I/O of files usually used to store radial velocity measurements.
    """
    def readuncommented(self):
        lines = self.readlines()
        uncommented = [l for l in lines if l.strip()[0].isdigit()]
        return uncommented