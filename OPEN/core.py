# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
from numpy import polyfit, RankWarning, append
# see http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
import warnings
from logger import clogger, logging
import sys

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


def do_restrict(system, quantity, maxval):
	if quantity == 'error':
		print 'the max is ', maxval

		# we have to keep a record of how many values come out of each file
		t, rv, err = system.time_full, system.vrad_full, system.error_full # temporaries
		for i, (fname, [n1, n2]) in enumerate(sorted(system.provenance.iteritems())):
			print i, n1, n2
			# print err[:n1]

			val = err[:n1] <= maxval
			nout = (val == False).sum()
			# system.provenance keeps the record
			if nout >= n1: 
				system.provenance[fname][1] = n1
			else:
				system.provenance[fname][1] = nout

			print (val == False).sum(), n1
			t, rv, err = t[n1:], rv[n1:], err[n1:]

		# now build the full boolean vector 
		vals = system.error_full <= maxval
		print vals, len(vals)

		# and pop out the values from time, vrad, and error
		# leaving all *_full vectors intact
		system.time = system.time_full[vals]
		system.vrad = system.vrad_full[vals]
		system.error = system.error_full[vals]

		# vals = system.error_full <= maxval
		# print 'From a total of ', len(system.error_full), \
		#       ' we are removing ', (vals == False).sum()

		# system.time = system.time_full[vals]
		# system.vrad = system.vrad_full[vals]
		# system.error = system.error_full[vals]

