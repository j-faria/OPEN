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
try:
    import triangle
except ImportError:
    import corner as triangle
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
    a, tau, gamma, period = p
    if (0.<a<1. and 5.<tau<50. and 20.<gamma<40. and 5.<period<50.): 
        return 0.0
    return -np.inf

def lnprob1(p, x, y, yerr):
    lp = lnprior1(p)
    return lp + lnlike1(p, x, y, yerr) if np.isfinite(lp) else -np.inf

GPfuncs['QuasiPeriodic'] = (k1, lnlike1, lnprior1, lnprob1)


k2 = lambda a,jit,tau,gamma,period: WhiteKernel(jit) + a * ExpSquaredKernel(tau) * ExpSine2Kernel(gamma, period)
kernel = k2(*[1., 1., 1., 1., 1.])
gp = george.GP(kernel, solver=george.HODLRSolver)
def lnlike2(p, t, y, yerr):
    # z1, z2, z3, z4 = p
    # kernel = z1**2 * kernels.ExpSquaredKernel(z2**2) * kernels.ExpSine2Kernel(2./z4**2, z3)
    
    # kernel = k2(*p)
    # gp = george.GP(kernel, solver=george.HODLRSolver)

    gp.kernel.pars = p
    gp.compute(t, yerr)
    # print gp.lnlikelihood(y)
    return gp.lnlikelihood(y)

def lnprior2(p):
    a, jit, tau, gamma, period = p
    if (0.0005 < a < 0.5 and 1e-7<jit<1e-4 and 40<tau<100000 and 0.2<gamma<5 and 10<period<40):
        return 0.0
    return -np.inf

def lnprob2(p, x, y, yerr):
    lp = lnprior2(p)
    return lp + lnlike2(p, x, y, yerr) if np.isfinite(lp) else -np.inf

GPfuncs['QuasiPeriodicJitter'] = (k2, lnlike2, lnprior2, lnprob2)


def fit_gp(model, initial, data, ncpu, nwalkers=20):
    k, lnlike, lnprior, lnprob = model
    ndim = len(initial)
    p0 = [np.array(initial) + rel(initial, 1) * np.random.randn(ndim) for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data, threads=1)

    msg = blue('    :: ') + 'Running burn-in...'
    clogger.info(msg)

    p0, lnp, _ = sampler.run_mcmc(p0, 100)
    sampler.reset()

    p0, lnp, _ = sampler.run_mcmc(p0, 200)
    sampler.reset()

    niter = 1000
    
    msg = blue('    :: ') + 'Running %d MCMC chains for %d iterations...' % (nwalkers, niter)
    clogger.info(msg)

    logl = []
    t1 = time()
    for p0, lnp, _ in sampler.sample(p0, None, None, iterations=niter):
    #     pass
        logl.append(max(lnp))
    # p0, lnp, _ = sampler.run_mcmc(p0, niter)
    t2 = time()
    logl = np.array(logl)

    p = p0[np.argmax(lnp)]


    msg = blue('    :: ') + 'MCMC took %f seconds' % (t2-t1)
    clogger.info(msg)
    return sampler, p, logl


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
    if training_variable == 'bis_span':
        yerr = 2.e-3*system.error

    
    # subtract mean
    y = y - np.mean(y)
    data = (t, y, yerr)

    model = GPfuncs['QuasiPeriodicJitter']

    # print y.ptp()
    initial = np.array([0.01, 1e-5, 5000, 1, 23])
    # best_p = initial
    sampler, best_p, logl = fit_gp(model, initial, data, ncpu)
    samples = sampler.flatchain 
    std = samples.std(axis=0)

    msg = yellow('    :: ') + 'Best GP hyperparameters: ' + initial.size*' %f ' % tuple(best_p)
    clogger.info(msg)
    msg = yellow('    :: ') + 'std of the chains:       ' + initial.size*' %f ' % tuple(std)
    clogger.info(msg)




    plt.figure()
    for i in range(samples.shape[1]+1):
        plt.subplot(6,1,i+1)
        if i == samples.shape[1]:
            plt.plot(logl)
        else:
            plt.plot(samples[:,i])
    plt.show()

    


    # # The positions where the prediction should be computed.
    x = np.linspace(min(t), max(t), 5000)
    x = np.hstack((x, t))
    x.sort()

    # # Plot 24 posterior samples.

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
    print gp.lnlikelihood(y)
    # Compute the prediction conditioned on the observations and plot it.
    # t1 = time()
    m, cov = gp.predict(y, x)
    m1, cov = gp.predict(y, t)
    # print time() - t1
    plt.figure()
    plt.subplot(211)

    # phase, fwhm_sim = np.loadtxt('/home/joao/phd/data/simulated/HD41248/HD41248_simul_oversampled.rdb', unpack=True, usecols=(0, 4), skiprows=2)
    # plt.plot(phase*18.3+t[0], fwhm_sim - fwhm_sim.mean(), 'g-')

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

    plt.figure()
    ax = plt.subplot(211)
    ts = BasicTimeSeries()
    ts.time = t
    ts.vrad = y
    ts.error = yerr
    per = gls(ts)
    per._plot(axes=ax, newFig=False)
    ax = plt.subplot(212)
    ts.vrad = y-m1
    per = gls(ts)
    per._plot(axes=ax, newFig=False)

    plt.show()
    # sys.exit(0)

    # fig = triangle.corner(samples, plot_contours=False)

    enter = raw_input('Press Enter to continue: ')
    if enter == 'n':
        sys.exit(0)

    return best_p, std
