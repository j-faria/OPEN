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
from numpy import loadtxt, savetxt, mean

# intra-package imports
from .logger import clogger, logging
from .utils import day2year, rms, ask_yes_no
from shell_colors import blue, yellow

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)


def read_rv(*filenames, **kwargs):
    """
    Read one or more files containing radial velocity measurements.
    
    Parameters
    ----------
    filenames: string
        One or more files to read.

    Optional kwargs
    ---------------
    skip: number of lines to skip in the files' headers
    format: column-format of the files
    verbose: verbosity toggle
    
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
    others: tuple
        Extra columns present in the file
    """

    # set logging level
    clogger.setLevel(logging.VERBOSE) \
                if (kwargs.has_key('verbose') and kwargs['verbose']) \
                else clogger.setLevel(logging.INFO)

    # how many header lines to skip?
    if (kwargs.has_key('skip')): header_skip = int(kwargs['skip'])
    # format of file
    if (kwargs.has_key('format')): format = kwargs['format']
    format = 'drs35' if (format is None) else format

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
    if format == 'drs35': # default
        t, rv, err, \
        fwhm, contrast, bis_span, noise, s_mw, sig_s, \
        rhk, sig_rhk, sn_CaII, sn10, sn50, sn60 = loadtxt(files, unpack=True)
        others = (fwhm, contrast, bis_span, noise, s_mw, sig_s, rhk, sig_rhk, sn_CaII, sn10, sn50, sn60)

    elif format == 'drs34' or format == 'coralie':
        t, rv, err, \
        fwhm, contrast, bis_span, noise, sn10, sn50, sn60  = loadtxt(files, unpack=True)
        others = (fwhm, contrast, bis_span, noise, sn10, sn50, sn60)

    # elif format == 'coralie':
    #     t, rv, err, 
    #     fwhm, contrast, bis_span, noise, sn10, sn50, sn60 = loadtxt(files, unpack=True)
    #     others = (fwhm, contrast, bis_span, noise, sn10, sn50, sn60)

    elif format == 'basic':
        t, rv, err = loadtxt(files, unpack=True, usecols=(0,1,2))
        others = ()

    
    # verbose stats about data
    info = blue('INFO: ') 
    sinfo = blue('    : ') 
    stats = None
    if (kwargs.has_key('verbose') and kwargs['verbose']):
        tspan = max(t) - min(t)
        rvspan = max(rv) - min(rv)
        stats = '\n'
        stats += info + "Timespan : %f days = %f years   ---   %fJD, %fJD\n" % (tspan, day2year(tspan), max(t), min(t))
        stats += sinfo + "RV span  : %f km/s = %f m/s\n" % (rvspan, rvspan*1e3)
        stats += sinfo + "RV rms [m/s] : %f\n\n" % rms(rv)
        stats += sinfo + "{:14s} : {:10.3f}\n".format('<RV> [km/s]', mean(rv))
        if format in ('drs35', 'drs34', 'coralie'):
            stats += sinfo + "{:14s} : {:10.3f}\n".format('<fwhm> [km/s]', mean(others[0]))
            stats += sinfo + "{:14s} : {:10.3f}\n".format('<contrast>', mean(others[1]))
            stats += sinfo + "{:14s} : {:10.3f}\n".format('<BIS> [km/s]', mean(others[2]))
            if format in ('drs35'):
                stats += sinfo + "{:14s} : {:10.3f}\n".format('<S_index> [MW]', mean(others[4]))
                stats += sinfo + "{:14s} : {:10.3f}\n".format('<log(rhk)>', mean(others[6]))
 
    clogger.verbose(stats)
    
    return t, rv, err, dic, others
    
    

def write_rv(system, filename, **kwargs):
    """
    Write system's RVs (and everything else) to a file.
    
    Parameters
    ----------
    system: instance of rvSeries
        System for which to save information.
    filename: string
        Name of the output file.

    Optional kwargs
    ---------------
    
    Returns
    -------
        nothing
    """
    print filename

    if os.path.isfile(filename):
        # this file exists
        clogger.warning(yellow('Warning: ')+'File already exists. Replace? [y/n]')
        if ask_yes_no(' '):
            pass
        else:
            clogger.info(blue('INFO: ')+'Aborting.')
            return

    header_basic1 = ['jdb', 'vrad', 'svrad']
    header_basic2 = ['-'*len(s) for s in header_basic1]

    header_drs35_1 = ['jdb', 'vrad', 'svrad', 'fwhm', 'contrast', 'bis_span', \
                     'noise', 's_mw', 'sig_s', 'rhk', 'sig_rhk', 'sn_CaII', \
                     'sn10', 'sn50', 'sn60']
    header_drs35_2 = ['-'*len(s) for s in header_drs35_1]

    header_drs34_1 = ['jdb', 'vrad', 'svrad', 'fwhm', 'contrast', 'bis_span', \
                      'noise', 'sn10', 'sn50', 'sn60']
    header_drs34_2 = ['-'*len(s) for s in header_drs34_1]

    if len(system.extras) == 0: # basic
        head = '\t'.join(header_basic1) + '\n' + '\t'.join(header_basic2)
        fmt = ['%12.6f', '%9.5f', '%7.5f']
        X = zip(system.time, system.vrad, system.error)
        savetxt('teste.rdb', X, fmt=fmt, \
                delimiter='\t', header=head, comments='')

    elif len(system.extras) == 7:
        head = '\t'.join(header_drs34_1) + '\n' + '\t'.join(header_drs34_2)
        fmt = ['%12.6f', '%9.5f', '%7.5f'] + ['%7.5f']*7
        X = zip(system.time, system.vrad, system.error, *list(system.extras))
        savetxt('teste.rdb', X, fmt=fmt, \
                delimiter='\t', header=head, comments='')

    elif len(system.extras) == 12:
        head = '\t'.join(header_drs35_1) + '\n' + '\t'.join(header_drs35_2)
        fmt = ['%12.6f', '%9.5f', '%7.5f'] + ['%7.5f']*12
        X = zip(system.time, system.vrad, system.error, *list(system.extras))
        savetxt('teste.rdb', X, fmt=fmt, \
                delimiter='\t', header=head, comments='')



class rvfile(file):
    """
    Subclass of Python's File class that implements specific methods
    for I/O of files usually used to store radial velocity measurements.
    """
    def readuncommented(self):
        lines = self.readlines()
        uncommented = [l for l in lines if l.strip()[0].isdigit()]
        return uncommented