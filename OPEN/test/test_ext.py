# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

def test_periodogram_ext():
	""" Test extensions related to periodograms """
	try:
		from ..ext import blombscargle
	except ImportError, e:
		raise e 