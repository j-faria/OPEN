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
