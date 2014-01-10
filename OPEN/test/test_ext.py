# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
from numpy import linspace, zeros_like, savetxt
import difflib


def test_periodogram_ext():
	""" Test extensions related to periodograms """
	try:
		from ..ext import blombscargle
	except ImportError, e:
		raise e 


def test_rv_ext():
	""" Test extensions related to building rv curves """
	# test imports
	try:
		from ..ext.get_rvN import get_rvn
		from ..ext.get_rv import get_rv
	except ImportError, e:
		raise e 

	times = linspace(2449460, 2452860, 100)
	vel1 = zeros_like(times)
	vel2 = zeros_like(times)

	get_rv(times, 1425, 10, 0.9, 0.2, 2452000, vel1)
	get_rvn(times, [1425, 13], [10, 3], [0.9, 0.02], [0.2, 0.3], [2452000, 2451000], 0., vel2)

	savetxt('out_test.txt', zip(times, vel1, vel2))

	f1 = open('out_normal.txt')
	f2 = open('out_test.txt')
	diff = difflib.SequenceMatcher(None, f1.read(), f2.read())
	assert diff.ratio() == 1., 'Non-normal RV curve building: something is wrong with get_rv or get_rvn extensions'