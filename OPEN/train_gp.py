# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

# standard library imports
import sys

# other imports
import matplotlib.pyplot as plt
import numpy as np
import george
from george import kernels
import triangle
import emcee
from time import time

# intra-package imports
from logger import clogger
from shell_colors import yellow, red, blue

def lnlike(p, t, y, yerr):
    z1, z2, z3, z4 = p
    kernel = z1**2 * kernels.ExpSquaredKernel(z2**2) * kernels.ExpSine2Kernel(2./z4**2, z3)
    gp = george.GP(kernel)
    gp.compute(t, yerr)
    return gp.lnlikelihood(y)

def lnprior(p):
    z1, z2, z3, z4 = p
    # if (0. < z1 < 0.01 and 5. < z2 < 50. and 20. < z3 < 30. and 5. < z4 < 50.): 
    if (0.001 < z1 < 0.009 and 20. < z2 < 30. and 20. < z3 < 30. and 25. < z4 < 40.): 
        return 0.0
    return -np.inf

def lnprob(p, x, y, yerr):
    lp = lnprior(p)
    return lp + lnlike(p, x, y, yerr) if np.isfinite(lp) else -np.inf


def fit_gp(initial, data, nwalkers=10):
    ndim = len(initial)
    p0 = [np.array(initial) #+ 1e-8 * np.random.randn(ndim)
          for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data)

    msg = blue('    :: ') + 'Running burn-in...'
    clogger.info(msg)

    p0, lnp, _ = sampler.run_mcmc(p0, 1000)
    sampler.reset()

    niter = 10
    
    msg = blue('    :: ') + 'Running %d MCMC chains for %d iterations...' % (nwalkers, niter)
    clogger.info(msg)

    t1 = time()
    p0, lnp, _ = sampler.run_mcmc(p0, niter)
    t2 = time()

    p = p0[np.argmax(lnp)]


    msg = blue('    :: ') + 'MCMC took %f seconds' % (t2-t1)
    clogger.info(msg)
    return sampler, p



def do_it(system, training_variable):

	t = system.time

	# find the quantity on which to train the GP
	i = system.extras._fields.index(training_variable) # index corresponding to the quantity
	y = system.extras[i]

	if training_variable == 'rhk': 
		training_variable_error = 'sig_rhk'
		i = system.extras._fields.index(training_variable_error) # index corresponding to the uncertainties
		yerr = system.extras[i]
	if training_variable == 'fwhm': 
		yerr = 2.35 * system.error
	
	
	# subtract mean
	y = y - np.mean(y)

	data = (t, y, 1.0 / yerr ** 2)
	initial = np.array([0.005, 27., 25.3, 30.])
	sampler, best_p = fit_gp(initial, data)

	msg = yellow('    :: ') + 'Best GP hyperparameters: z1=%f, z2=%f, z3=%f, z4=%f' % tuple(best_p)
	clogger.info(msg)


	


	# # The positions where the prediction should be computed.
	x = np.linspace(min(t), max(t), 1000)
	x = np.hstack((x, t))
	x.sort()

	# # Plot 24 posterior samples.
	# # samples = sampler.flatchain
	# # for s in samples[np.random.randint(len(samples), size=4)]:
	# #     # Set up the GP for this sample.
	# #     z1, z2, z3, z4 = s
	# #     kernel = z1**2 * kernels.ExpSquaredKernel(z2**2) * kernels.ExpSine2Kernel(2./z4**2, z3)
	# #     gp = george.GP(kernel)
	# #     gp.compute(t, yerr)

	# #     # Compute the prediction conditioned on the observations and plot it.
	# #     m = gp.sample_conditional(y, x)
	# #     plt.plot(x, m, color="#4682b4", alpha=0.3)

	# plot lnp solution
	z1, z2, z3, z4 = best_p
	kernel = z1**2 * kernels.ExpSquaredKernel(z2**2) * kernels.ExpSine2Kernel(2./z4**2, z3)
	gp = george.GP(kernel)
	gp.compute(t, yerr)
	# Compute the prediction conditioned on the observations and plot it.
	# t1 = time()
	m = gp.sample_conditional(y, x)
	m1 = gp.sample_conditional(y, t)
	# print time() - t1
	plt.figure()
	plt.subplot(211)
	plt.plot(x, m, color='r', alpha=0.8)
	plt.plot(t, m1, color='r', alpha=0.8)

	# Plot the data
	plt.errorbar(t, y, yerr=yerr, fmt=".k", capsize=0)
	plt.ylabel(training_variable)

	# Plot the residuals
	plt.subplot(212)
	plt.errorbar(t, y - m1, yerr=yerr, fmt=".k", capsize=0)


	plt.xlabel('Time [days]')

	plt.show()

	# sys.exit(0)
	return best_p
	# fig = triangle.corner(samples)