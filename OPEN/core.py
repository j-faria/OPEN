# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
from numpy import polyfit, RankWarning
# see http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
import warnings
from logger import clogger, logging

def do_fit(system, verbose):
	from shell_colors import yellow
	degree = system.model['d']
	with warnings.catch_warnings(record=True) as w:
		p = polyfit(system.time, system.vrad, degree)
		if len(w):
			msg = yellow('Warning: ') + 'Polyfit may be poorly conditioned. ' \
			      + 'Maybe try a lower degree drift?'
			clogger.info(msg)
	return p