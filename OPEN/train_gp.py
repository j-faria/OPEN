# -*- coding: utf-8 -*-
#
# Copyright (c) 2013-2015 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

"""
This module implements Gaussian Process regression on indicators like FWHM or logRhk.
After the covariance structure in one of these variables is determined, it can be used
for the radial velocities as well, by considering a GP with a keplerian mean function.
"""
# standard library imports
import sys

# other imports
import matplotlib.pyplot as plt
import numpy as np
import george
from george import kernels
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel
import triangle
import emcee
from time import time

# intra-package imports
from logger import clogger
from shell_colors import yellow, red, blue
from classes import BasicTimeSeries
from periodograms import gls


def rel(array, percentage):
    """ Get relative value from array """
    return (percentage/100.) * array




# We implement different GP covariances
GPfuncs = {}

k1 = lambda a,tau,gamma,period: a**2 * ExpSquaredKernel(tau**2) * ExpSine2Kernel(2./gamma**2, period)
def lnlike1(p, t, y, yerr):
    # z1, z2, z3, z4 = p
    # kernel = z1**2 * kernels.ExpSquaredKernel(z2**2) * kernels.ExpSine2Kernel(2./z4**2, z3)
    kernel = k1(*p)
    gp = george.GP(kernel, solver=george.HODLRSolver)
    gp.compute(t, yerr)
    return gp.lnlikelihood(y)

def lnprior1(p):
    z1, z2, z3, z4 = p
    if (0.<z1<1. and 5.<z2<50. and 20.<z3<40. and 5.<z4<50.): 
        return 0.0
    return -np.inf

def lnprob1(p, x, y, yerr):
    lp = lnprior1(p)
    return lp + lnlike1(p, x, y, yerr) if np.isfinite(lp) else -np.inf

GPfuncs['QuasiPeriodic'] = (k1, lnlike1, lnprior1, lnprob1)


k2 = lambda a,jit,tau,gamma,period: WhiteKernel(jit) + a * ExpSquaredKernel(tau) * ExpSine2Kernel(gamma, period)
def lnlike2(p, t, y, yerr):
    # z1, z2, z3, z4 = p
    # kernel = z1**2 * kernels.ExpSquaredKernel(z2**2) * kernels.ExpSine2Kernel(2./z4**2, z3)
    kernel = k2(*p)
    gp = george.GP(kernel, solver=george.HODLRSolver)
    gp.compute(t, yerr)
    return gp.lnlikelihood(y)

def lnprior2(p):
    a, jit, tau, gamma, period = p
    if (0.0005 < a < 0.005 and 1e-6<jit<1e-4 and 50<tau<90 and 0.2<gamma<5 and 25<period<35):
        return 0.0
    return -np.inf

def lnprob2(p, x, y, yerr):
    lp = lnprior2(p)
    return lp + lnlike2(p, x, y, yerr) if np.isfinite(lp) else -np.inf

GPfuncs['QuasiPeriodicJitter'] = (k2, lnlike2, lnprior2, lnprob2)


def fit_gp(model, initial, data, ncpu, nwalkers=10):
    k, lnlike, lnprior, lnprob = model
    ndim = len(initial)
    p0 = [np.array(initial) + rel(initial, 1) * np.random.randn(ndim) for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data, threads=1)

    msg = blue('    :: ') + 'Running burn-in...'
    clogger.info(msg)

    p0, lnp, _ = sampler.run_mcmc(p0, 100)
    sampler.reset()

    niter = 250
    
    msg = blue('    :: ') + 'Running %d MCMC chains for %d iterations...' % (nwalkers, niter)
    clogger.info(msg)

    t1 = time()
    p0, lnp, _ = sampler.run_mcmc(p0, niter)
    t2 = time()

    p = p0[np.argmax(lnp)]


    msg = blue('    :: ') + 'MCMC took %f seconds' % (t2-t1)
    clogger.info(msg)
    return sampler, p


def do_it(system, training_variable, ncpu=1):

    t = system.time

    # find the quantity on which to train the GP
    i = system.extras._fields.index(training_variable) # index corresponding to the quantity
    y = system.extras[i]

    if training_variable == 'rhk': 
        training_variable_error = 'sig_rhk'
        i = system.extras._fields.index(training_variable_error) # index corresponding to the uncertainties
        yerr = system.extras[i]
    if training_variable == 'fwhm':
        if system.units == 'm/s':
            f = 2.35e-3
        else:
            f = 2.35
        yerr = f * system.error
    
    # subtract mean
    y = y - np.mean(y)
    data = (t, y, 1.0 / yerr ** 2)

    model = GPfuncs['QuasiPeriodicJitter']
    initial = np.array([0.001, 1e-5, 80, 1, 30])
    sampler, best_p = fit_gp(model, initial, data, ncpu)

    msg = yellow('    :: ') + 'Best GP hyperparameters: ' + initial.size*' %f ' % tuple(best_p)
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
    best_p[1] = 0.
    kernel = model[0](*best_p)
    gp = george.GP(kernel, solver=george.HODLRSolver)
    gp.compute(t, yerr)
    # Compute the prediction conditioned on the observations and plot it.
    # t1 = time()
    m = gp.sample_conditional(y, x)
    m1 = gp.sample_conditional(y, t)
    # print time() - t1
    plt.figure()
    plt.subplot(211)
    plt.plot(x, m, color='r', alpha=0.8)
    # plt.plot(t, m1, color='r', alpha=0.8)
    # Plot the data
    plt.errorbar(t, y, yerr=yerr, fmt=".k", capsize=0)
    # plt.plot(t, system.extras.rhk_activity - system.extras.rhk_activity.mean(), "og")
    plt.ylabel(training_variable)
    # Plot the residuals
    plt.subplot(212)
    plt.errorbar(t, y - m1, yerr=yerr, fmt=".k", capsize=0)

    plt.xlabel('Time [days]')

    ts = BasicTimeSeries()
    ts.time = t
    ts.vrad = y-m1
    ts.error = yerr
    per = gls(ts)
    per._plot()

    # plt.show()
    # sys.exit(0)
    enter = raw_input('Press Enter to continue: ')
    if enter == 'n':
        sys.exit(0)

    return best_p
    # fig = triangle.corner(samples)