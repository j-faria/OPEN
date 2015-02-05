# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

"""
This module contains the analytical expressions for the priors used in the
MCMC/MultiNest analyses in OPEN.
"""

import numpy as np
from scipy.stats import beta as bt

def uniform(x, u1, u2):
	""" Normalized uniform prior """
	range_uniform = u2-u1
	return np.ones_like(x) / range_uniform


def jeffreys(x, jmin, jmax):
	""" Normalized jeffreys prior """
	f = np.log(jmax/jmin)
	return 1./(x*f)


def modjeffreys(x, mj0, mjmax):
	""" Normalized modified jeffreys prior """
	f = np.log(1 + mjmax/mj0)
	return 1./((x+mj0)*f)


def beta(x, a, b):
	""" Beta prior for shape parameters a and b """
	return bt.pdf(x, a, b)
