# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

# standard library imports
import warnings
import datetime, time
import subprocess
import sys, os
import random
from collections import namedtuple
from string import ascii_lowercase
import copy
import fileinput
import tempfile

# other imports
from numpy import polyfit, RankWarning, append, zeros_like, savetxt
import numpy as np
import pylab
import matplotlib.pyplot as plt
from deap import base, creator, tools, algorithms

# see http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
from scipy.optimize import leastsq
from scipy.stats.stats import spearmanr, pearsonr
from scipy.stats import nanmean, nanstd

# intra-package imports
from .classes import MCMC
from logger import clogger, logging
from ext.get_rvN import get_rvn
from .periodograms import ls_PressRybicki, gls
try:
	from ext.periodogram_DF import periodogram_DF
	periodogram_DF_available = True
except ImportError:
	periodogram_DF_available = False
from shell_colors import yellow, red, blue
from .utils import julian_day_to_date, ask_yes_no

timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
pi = np.pi

def do_fit(system, verbose):
	try:
		degree = system.model['d']
	except TypeError:
		msg = yellow('Warning: ') + 'To remove linear trends, run mod before fit.\n'
		clogger.error(msg)

	# with warnings.catch_warnings(record=True) as w:
	# 	p = polyfit(system.time, system.vrad, degree)
	# 	if len(w):
	# 		msg = yellow('Warning: ') + 'Polyfit may be poorly conditioned. ' \
	# 		      + 'Maybe try a lower degree drift?'
	# 		clogger.info(msg)
	# return p
	try:
		system.per
	except AttributeError:
		system.per = ls_PressRybicki(system)  # this periodogram is very fast
	finally:
		peak = system.per.get_peaks(output_period=True)

	peak_val = system.per.power.max()
	peak_fap = system.per.fap
	if peak_fap < 0.01:
		msg = blue('INFO: ') + 'Peak found at %f days with FAP ~ %e' % (peak[1], peak_fap)
		clogger.info(msg)
	else:
		msg = yellow('Warning: ') + 'No peaks found with FAP < 1%'
		clogger.error(msg)
		return None

	p0 = peak[1] # initial guess period
	k0 = system.vrad.max()
	e0 = 0.2
	om0 = 10
	t0 = system.time.min()
	g0 = system.vrad.mean()
	msg = blue('INFO: ') + 'Setting initial guesses' #, p0, k0, e0, om0, t0, g0
	initial_guess = [p0, k0, e0, om0, t0, g0]
	## call levenberg markardt fit
	lm = do_lm(system, initial_guess)
	lm_par = lm[0]
	
	## loop over planets
	msg = yellow('RESULT: ') + 'Best fit is'
	clogger.info(msg)
	print("%3s %12s %10s %10s %10s %15s %9s" % \
		('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
	for i, planet in enumerate(list(ascii_lowercase)[:1]):
		P, K, ecc, omega, T0, gam = [lm_par[j::6] for j in range(6)]
		print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )
	
	# save fit in the system
	system.save_fit(lm_par, 0.)
	par1 = lm_par # backup

	j=1
	while j<4:
		# is there a relevant peak in the periodogram of the residuals?
		system.per_resid = gls(system, quantity='resid')
		# system.per_resid._plot()
		peak = system.per_resid.get_peaks(output_period=True)
		peak_val = system.per_resid.power.max()
		peak_fap = system.per_resid.FAP(peak_val)
		if peak_fap < 0.01:
			msg = blue('INFO: ') + 'Peak found at %f days with FAP ~ %e' % (peak[1], peak_fap)
			clogger.info(msg)
		else:
			msg = yellow('Warning: ') + 'No peaks found in the periodogram of the residuals with FAP < 1%'
			clogger.error(msg)
			return None

		new_system = copy.deepcopy(system)
		new_system.vrad = new_system.fit['residuals']

		p0 = peak[1] # initial guess period
		k0 = new_system.vrad.max()
		e0 = 0.2
		om0 = 10
		t0 = new_system.time.min()
		g0 = new_system.vrad.mean()
		msg = blue('INFO: ') + 'Setting initial guesses' #, p0, k0, e0, om0, t0, g0
		initial_guess = [p0, k0, e0, om0, t0, g0]
		## call levenberg markardt fit
		lm = do_lm(new_system, initial_guess)
		lm_par = lm[0]
		
		## loop over planets
		msg = yellow('RESULT: ') + 'Best fit is'
		clogger.info(msg)
		print("%3s %12s %10s %10s %10s %15s %9s" % \
			('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
		for i, planet in enumerate(list(ascii_lowercase)[:1]):
			P, K, ecc, omega, T0, gam = [lm_par[j::6] for j in range(6)]
			print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )

		par_all = np.append(par1, lm_par)
		# save fit in the system
		system.save_fit(par_all, 0.)
		j+=1

	return None


def do_restrict(system, quantity, *args):
	## restrict by uncertainty value
	if quantity == 'error':
		msg = blue('INFO: ') + 'Removing data with uncertainty higher than %f km/s' % args[0]
		clogger.info(msg)
		maxerr = args[0]
		return

		# we have to keep a record of how many values come out of each file
		t, rv, err = system.time_full, system.vrad_full, system.error_full # temporaries
		for i, (fname, [n1, n2]) in enumerate(sorted(system.provenance.iteritems())):
			print i, n1, n2
			# print err[:n1]

			val = err[:n1] <= maxerr
			nout = (val == False).sum()
			# system.provenance keeps the record
			if nout >= n1: 
				system.provenance[fname][1] = n1
			else:
				system.provenance[fname][1] = nout

			print (val == False).sum(), n1
			t, rv, err = t[n1:], rv[n1:], err[n1:]

		# now build the full boolean vector 
		vals = system.error_full <= maxerr
		print vals, len(vals)

		# and pop out the values from time, vrad, and error
		# leaving all *_full vectors intact
		system.time = system.time_full[vals]
		system.vrad = system.vrad_full[vals]
		system.error = system.error_full[vals]


	## restrict by date (JD)
	if quantity == 'date':
		msg = blue('INFO: ') + 'Retaining data between %i and %i JD' % (args[0], args[1])
		clogger.info(msg)
		minjd, maxjd = args[0], args[1]

		# we have to keep a record of how many values come out of each file
		t, rv, err = system.time_full, system.vrad_full, system.error_full # temporaries
		for i, (fname, [n1, n2]) in enumerate(sorted(system.provenance.iteritems())):
			lower = t[:n1] <= minjd
			higher = t[:n1] >= maxjd
			keep = ~(lower | higher)
			nout = (keep == False).sum()
			# system.provenance keeps the record
			if nout >= n1: 
				system.provenance[fname][1] = n1
			else:
				system.provenance[fname][1] = nout

			t, rv, err = t[n1:], rv[n1:], err[n1:]

		# now build the full boolean vector 
		lowers = system.time_full <= minjd
		highers = system.time_full >= maxjd
		# fancy syntax just to negate the intersection of lowers and highers
		keepers = ~(lowers | highers)
		print keepers

		# and pop out the values from time, vrad, and error
		# leaving all *_full vectors intact
		system.time = system.time_full[keepers]
		system.vrad = system.vrad_full[keepers]
		system.error = system.error_full[keepers]	
		# also from extras, but this is trickier
		d = system.extras._asdict() # asdict because namedtuple is immutable
		for i, field in enumerate(system.extras._fields):
			d[field] = system.extras_full[i][keepers]
		extra = namedtuple('Extra', system.extras_names, verbose=False)
		system.extras = extra(**d)

	## restrict to values from one year
	if quantity == 'year':
		msg = blue('INFO: ') + 'Retaining data from %i' % args[0]
		clogger.info(msg)
		yr = args[0]

		# we have to keep a record of how many values come out of each file
		time = system.time_full # temporaries
		for i, (fname, [n1, n2]) in enumerate(sorted(system.provenance.iteritems())):
			years = np.array([julian_day_to_date(t)[0] for t in time[:n1]])
			keepers = years == yr
			nout = (keepers == True).sum()
			# system.provenance keeps the record
			if nout >= n1: 
				system.provenance[fname][1] = n1
			else:
				system.provenance[fname][1] = n1 - nout
			time = time[n1:]

		# build the full boolean vector
		years = np.array([julian_day_to_date(t)[0] for t in system.time_full])
		keepers = years == yr

		# and pop out the values from time, vrad, and error
		# leaving all *_full vectors intact
		system.time = system.time_full[keepers]
		system.vrad = system.vrad_full[keepers]
		system.error = system.error_full[keepers]
		# also from extras
		d = system.extras._asdict() # asdict because namedtuple is immutable
		for i, field in enumerate(system.extras._fields):
			d[field] = system.extras_full[i][keepers]
		extra = namedtuple('Extra', system.extras_names, verbose=False)
		system.extras = extra(**d)	

	## restrict to values from a year range
	if quantity == 'years':
		msg = blue('INFO: ') + 'Retaining data between %i and %i' % (args[0], args[1])
		clogger.info(msg)
		yr1, yr2 = args[0], args[1]

		# we have to keep a record of how many values come out of each file
		time = system.time_full # temporaries
		for i, (fname, [n1, n2]) in enumerate(sorted(system.provenance.iteritems())):
			years = np.array([julian_day_to_date(t)[0] for t in time[:n1]])
			keepers = (years >= yr1) & (years <= yr2)
			nout = (keepers == True).sum()
			# system.provenance keeps the record
			if nout >= n1: 
				system.provenance[fname][1] = n1
			else:
				system.provenance[fname][1] = n1 - nout
			time = time[n1:]


		# now build the full boolean vector 
		years = np.array([julian_day_to_date(t)[0] for t in system.time_full])
		keepers = (years >= yr1) & (years <= yr2)

		# and pop out the values from time, vrad, and error
		# leaving all *_full vectors intact
		system.time = system.time_full[keepers]
		system.vrad = system.vrad_full[keepers]
		system.error = system.error_full[keepers]	
		# also from extras, but this is trickier
		d = system.extras._asdict() # asdict because namedtuple is immutable
		for i, field in enumerate(system.extras._fields):
			d[field] = system.extras_full[i][keepers]
		extra = namedtuple('Extra', system.extras_names, verbose=False)
		system.extras = extra(**d)

	return

def do_demc(system, burnin=500):
	import demc.dream as dream
	from time import sleep
	from re import findall

	try:
		degree = system.model['d']
		keplerians = system.model['k']
	except TypeError:
		msg = red('Error: ') + 'Need to run mod before de. '
		clogger.error(msg)
		return	

	msg = blue('INFO: ') + 'Transfering data to DREAM...'
	clogger.info(msg)
	msg = blue('    : ') + 'Model is: %d keplerians + %d drift' % (keplerians, degree)
	clogger.info(msg)

	## write data to file to be read by DREAM
	dream_filename = 'input.rv'
	dream_header = 'file automatically generated for OPEN-DREAM analysis, ' + timestamp
	dream_header += '\n' + str(len(system.time))
	savetxt(dream_filename, zip(system.time, system.vrad, system.error),
		    header=dream_header,
		    fmt=['%12.6f', '%7.5f', '%7.5f'])


	## automatically fill the necessary parameters in the inlist
	replacer = '\tpar_num = %d\n' % (5*keplerians+1)
	for line in fileinput.input('OPEN/demc/namelist1', inplace=True):
		if 'par_num' in line: 
			print replacer,
		else:
			 print line,

	## get number of chains and generations from namelist
	with open('OPEN/demc/namelist1') as f:
		l1 = [line for line in f.readlines() if 'gen_num' in line][0]
		ngen = int(''.join(findall(r'[0-9]', l1)))
		f.seek(0)
		l2 = [line for line in f.readlines() if 'chain_num' in line][0]
		nc = int(''.join(findall(r'[0-9]', l2)))

	msg = blue('    : ') + 'Going to evolve %d chains for %d generations...' % (nc, ngen)
	clogger.info(msg)

	# let the user read the previous messages :)
	sleep(1.5)

	msg = blue('    : ') + 'Starting DREAM...'
	clogger.info(msg)

	dream.main()

	msg = blue('INFO: ') + 'DREAM: normal end of execution'
	clogger.info(msg)
	msg = blue('    : ') + 'Starting analysis of output...'
	clogger.info(msg)

	results = MCMC('OPEN/demc/problem1_chain*', burnin=burnin)
	results.read_chains()

	## output best solution
	best_ind = np.argmax(results.chains[1, :])
	par_best = results.chains[2:, best_ind]
	print # newline
	msg = yellow('RESULT: ') + 'Best solution is'
	clogger.info(msg)
	## loop over planets
	print("%3s %12s %10s %10s %10s %15s %9s" % \
		('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
	for i, planet in enumerate(list(ascii_lowercase)[:keplerians]):
		P, K, ecc, omega, T0, gam = [par_best[j::6] for j in range(6)]
		print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )
	
	print # newline
	msg = blue('INFO: ') + 'Estimating kernels...'
	clogger.info(msg)

	# p = results.kde1d('period_0', npoints=300)
	# ind = np.linspace(0.2, 5000, 300)

	results.do_plot_some_solutions(system)

	return results


def do_diffevol(system, just_de=False, npop=100, ngen=250):
	""" Carry out the fit using a (experimental) differential evolution
	algorithm """
	from de import DiffEvol

	try:
		degree = system.model['d']
		keplerians = system.model['k']
	except TypeError:
		msg = red('Error: ') + 'Need to run mod before de. '
		clogger.error(msg)
		return	

	# should we do it this way?
	system.vrad = system.vrad_full

	msg = blue('INFO: ') + 'Initializing DE algorithm...'
	clogger.info(msg)
	msg = blue('    : ') + 'Model is: %d keplerians + %d drift' % (keplerians, degree)
	clogger.info(msg)

	vel = zeros_like(system.time)

	## trend removal
	if degree > 0:
		z = np.polyfit(system.time, system.vrad, degree)
		poly = np.poly1d(z)
		system.vrad = system.vrad - poly(system.time)
		# save info to system
		system.model['drift'] = z
		if keplerians == 0: 
			## output some information
			msg = yellow('RESULT: ') + 'Fit of degree %d done! Coefficients:' % degree
			clogger.info(msg)
			msg = yellow('      : ') + str(z)
			clogger.info(msg)
			return

	## fitness function - N planet model
	def chi2_n(individual):
		""" Fitness function for N planet model """
		P, K, ecc, omega, T0, gam = [individual[i::6] for i in range(6)]
		#print ecc
		get_rvn(system.time, P, K, ecc, omega, T0, gam[0], vel)
		#print 'out of get_rvn'
		chi2 = sum(((system.vrad - vel)/system.error)**2)
		#print chi2
		return chi2

	pars = ['P', 'K', 'ecc', 'w', 'T0', 'gama']
	## parameter bounds 
	par_bounds = []
	par_bounds.append([5, 1000]) # period
	par_bounds.append([0.5, 150]) # semi amplitude
	par_bounds.append([0, 0.9]) # eccentricity
	par_bounds.append([0, 360]) # omega
	par_bounds.append([2350000, 2550000]) # periastron passage
	par_bounds.append([-100, 100]) # offset

	par_bounds = keplerians*par_bounds  # same bounds for each planet

	npar = 5*keplerians+1
	degfree = len(system.time) - npar
	
	de = DiffEvol(chi2_n, par_bounds, npop, C=0.1)

	msg = blue('INFO: ') + 'Created population with N=%d. Going to evolve for %d generations...' % (npop,ngen)
	clogger.info(msg)

	# res = de.optimize(ngen=ngen)
	print # newline
	for niter, res in enumerate(de(ngen=ngen)):
		if niter % 50 == 0:
			for i in range(keplerians):
				print('  ' + '\t'.join('{}: {:8.2f}'.format(*k) for k in zip(pars, de.minimum_location[6*i:6*(i+1)])))
			print '  Best chi^2: %5.2f \t chi2_red: %5.2f \n' % (de.minimum_value, de.minimum_value/degfree)

	best_param = de.minimum_location
	best_fitness = de.minimum_value

	## output results information
	msg = yellow('RESULT: ') + 'Best individual is'
	clogger.info(msg)

	## loop over planets
	print("%3s %12s %10s %10s %10s %15s %9s" % \
		('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
	for i, planet in enumerate(list(ascii_lowercase)[:keplerians]):
		P, K, ecc, omega, T0, gam = [best_param[j::6] for j in range(6)]
		# print P, K, ecc, omega, T0, gam
		print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )
	print #newline
	
	msg = yellow('RESULT: ') + 'Best fitness value: %7.2f  (%6.2f)\n' % (best_fitness, best_fitness/degfree)
	clogger.info(msg)

	if just_de:
		# save fit in the system
		system.save_fit(best_param, best_fitness)
		return

	msg = blue('INFO: ') + 'Calling LM to improve result...'
	clogger.info(msg)

	## call levenberg markardt fit
	lm = do_lm(system, best_param)
	lm_par = lm[0]
	
	## loop over planets
	msg = yellow('RESULT: ') + 'Best fit is'
	clogger.info(msg)
	print("%3s %12s %10s %10s %10s %15s %9s" % \
		('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
	for i, planet in enumerate(list(ascii_lowercase)[:keplerians]):
		P, K, ecc, omega, T0, gam = [lm_par[j::6] for j in range(6)]
		print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )
	print #newline

	chi2 = chi2_n(lm_par)
	msg = yellow('RESULT: ') + 'Best fitness value: %7.2f  (%6.2f)\n' % (chi2, chi2/(len(system.time)-npar))
	clogger.info(msg)

	# save fit in the system
	system.save_fit(lm_par, chi2)

def do_genetic(system, just_gen=False):
	""" Carry out the fit using a genetic algorithm and if 
	just_gen=False try to improve it with a run of the LM algorithm """
	try:
		degree = system.model['d']
		keplerians = system.model['k']
	except TypeError:
		msg = red('Error: ') + 'Need to run mod before gen. '
		clogger.error(msg)
		return

	# maxP = system.per.get_peaks(output_period=True)[1]
	# size_maxP = 10**(len(str(int(maxP)))-1)

	msg = blue('INFO: ') + 'Initializing genetic algorithm...'
	clogger.info(msg)
	msg = blue('    : ') + 'Model is: %d keplerians + %d drift' % (keplerians, degree)
	clogger.info(msg)

	vel = zeros_like(system.time)

	def chi2_1(individual):
		""" Fitness function for 1 planet model """
		P, K, ecc, omega, T0, gam = individual 
		get_rvn(system.time, P, K, ecc, omega, T0, gam, vel)
		chi2 = sum(((system.vrad - vel)/system.error)**2)
		#print chi2
		return chi2,

	def chi2_n(individual):
		""" Fitness function for N planet model """
		P, K, ecc, omega, T0, gam = [individual[i::6] for i in range(6)]
		#print ecc
		get_rvn(system.time, P, K, ecc, omega, T0, gam[0], vel)
		#print 'out of get_rvn'
		chi2 = sum(((system.vrad - vel)/system.error)**2)
		#print chi2
		return chi2,

	## create the required types -- the fitness and the individual.
	creator.create("FitnessMin", base.Fitness, weights=(-1.0,)) # minimization of a single objective
	creator.create("Individual", list, fitness=creator.FitnessMin) 

	## create parameters by sampling from their priors
	def P_prior():
		return random.uniform(5, 1000)
		# return random.gauss(maxP, size_maxP)
	def K_prior():
		return random.uniform(0, 150)
	def ecc_prior():
		return random.uniform(0, 0.9)
	def om_prior():
		return random.uniform(0, 360)
	def t0_prior():
		return random.uniform(2350000, 2550000)
	def gamma_prior():
		return random.uniform(-100, 100)
	priors = [P_prior, K_prior, ecc_prior, om_prior, t0_prior, gamma_prior]

	toolbox = base.Toolbox()
	toolbox.register("individual", tools.initCycle, creator.Individual, priors, n=keplerians)
	toolbox.register("population", tools.initRepeat, list, toolbox.individual)

	def mutPrior(individual, indpb):
		for i, fcn in enumerate(zip(individual, priors)):
			if random.random() < indpb:			
				individual[i] = fcn[1]()
		return individual,

	toolbox.register("evaluate", chi2_n)
	toolbox.register("mate", tools.cxTwoPoints)
	toolbox.register("mutate", mutPrior, indpb=0.10)
	toolbox.register("select", tools.selTournament, tournsize=3)

	npop = 500
	ngen = 150
	npar = 5*keplerians+1
	## build the population
	pop = toolbox.population(n=npop)
	## helper functions
	hof = tools.HallOfFame(1)
	stats = tools.Statistics(lambda ind: ind.fitness.values)
	stats.register("avg", nanmean)
	# stats.register("std", nanstd)
	stats.register("min", np.nanmin)
	# stats.register("max", np.nanmax)
	# stats.register("total", sigma3)
	stats.register("red", lambda v: min(v)/(len(system.time)-npar) )

	msg = blue('INFO: ') + 'Created population with N=%d. Going to evolve for %d generations...' % (npop,ngen)
	clogger.info(msg)

	## run the genetic algorithm
	algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=ngen, stats=stats, halloffame=hof, verbose=True)
	
	## output results information
	msg = yellow('RESULT: ') + 'Best individual is'
	clogger.info(msg)
	## loop over planets
	print("%3s %12s %10s %10s %10s %15s %9s" % \
		('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
	for i, planet in enumerate(list(ascii_lowercase)[:keplerians]):
		P, K, ecc, omega, T0, gam = [hof[0][j::6] for j in range(6)]
		print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )
	
	msg = yellow('RESULT: ') + 'Best fitness value: %s\n' % (hof[0].fitness)
	clogger.info(msg)

	if just_gen: 
		# save fit in the system and return, no need for LM
		system.fit['params'] = hof[0]
		system.fit['chi2'] = hof[0].fitness/(len(system.time)-npar)
		return

	msg = blue('INFO: ') + 'Calling LM to improve result...'
	clogger.info(msg)	

	## call levenberg markardt fit
	lm = do_lm(system, [hof[0][j::6] for j in range(6)])
	lm_par = lm[0]
	
	## loop over planets
	msg = yellow('RESULT: ') + 'Best fit is'
	clogger.info(msg)
	print("%3s %12s %10s %10s %10s %15s %9s" % \
		('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
	for i, planet in enumerate(list(ascii_lowercase)[:keplerians]):
		P, K, ecc, omega, T0, gam = [lm_par[j::6] for j in range(6)]
		print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )

	chi2 = chi2_n(lm_par)[0]
	msg = yellow('RESULT: ') + 'Best fitness value: %f, %f' % (chi2, chi2/(len(system.time)-npar))
	clogger.info(msg)

	# save fit in the system
	system.save_fit(lm_par, chi2)
	# system.fit['params'] = lm_par
	# system.fit['chi2'] = chi2

	# #  print p.minFitness, p.maxFitness, p.avgFitness, p.sumFitness
	# print 'Genetic:', p.bestFitIndividual, p.bestFitIndividual.fitness
	# lm = do_lm(system, p.bestFitIndividual.genes)
	# lm_par = lm[0]
	# print 'LM:', lm_par

	# # get best solution curve
	# new_time = system.get_time_to_plot()
	# vel = zeros_like(new_time)

	# P, K, ecc, omega, t0 = p.bestFitIndividual.genes
	# get_rv(new_time, P, K, ecc, omega, t0, vel)
	# # plot RV with time
	# plot(system.time, system.vrad, 'o')
	# plot(new_time, vel, '-')

	# P, K, ecc, omega, t0 = lm_par
	# get_rv(new_time, P, K, ecc, omega, t0, vel)
	# plot(new_time, vel, 'r-')
	# show()

	return

def do_lm(system, x0):

	vel = zeros_like(system.time)
	# print x0
	# print np.transpose(x0)

	def chi2_n_leastsq(params):
		""" Fitness function for N planet model """
		#print params
		P, K, ecc, omega, T0, gam = [params[i::6] for i in range(6)]
		#print ecc
		if any(e>1 or e<0 for e in ecc): return 1e99
		get_rvn(system.time, P, K, ecc, omega, T0, gam[0], vel)
		#print 'out of get_rvn'
		return (system.vrad - vel) / system.error

	def chi2_2(params):
		P, K, ecc, omega, T0, gam = params
		get_rvn(system.time, P, K, ecc, omega, T0, gam, vel)
		return system.vrad - vel

	# x0 = np.transpose(x0)
	x0 = np.abs(x0)
	return leastsq(chi2_n_leastsq, x0, full_output=0, ftol=1e-15)#, maxfev=10)


def do_multinest(system):
	msg = blue('INFO: ') + 'Transfering data to MultiNest...'
	clogger.info(msg)

	# write data to file to be read by MultiNest
	nest_filename = 'input.rv'
	nest_header = 'file automatically generated for MultiNest analysis, ' + timestamp
	nest_header += '\n' + str(len(system.time))
	savetxt(nest_filename, zip(system.time, system.vrad, system.error),
		    header=nest_header,
		    fmt=['%12.6f', '%7.5f', '%7.5f'])

	msg = blue('INFO: ') + 'Starting MultiNest...'
	clogger.info(msg)
	
	cmd = 'mpirun -np 2 ./OPEN/multinest/nest'
	subprocess.call(cmd, shell=True)
	# os.system(cmd)

	return

def do_correlate(system, vars=(), verbose=False):
	# just to be sure, but this should not pass through docopt in commands.py
	if len(vars) != 2: return

	var1 = vars[0]
	var2 = vars[1]

	# handle inexistent fields
	available = system.extras._fields + ('vrad',)
	if var1 not in available:
		msg = red('ERROR: ') + 'The name "%s" is not available for correlation.\n' % var1
		clogger.fatal(msg)
		return
	if var2 not in available:
		msg = red('ERROR: ') + 'The name "%s" is not available for correlation.\n' % var2
		clogger.fatal(msg)
		return

	if var1 == 'vrad': 
		v1 = system.vrad
	else:
		i = system.extras._fields.index(var1)
		v1 = system.extras[i]

	if var2 == 'vrad': 
		v2 = system.vrad
	else:
		i = system.extras._fields.index(var2)
		v2 = system.extras[i]		

	pr = pearsonr(v1, v2)
	sr = spearmanr(v1, v2)

	if verbose:
		print blue('[Pearson correlation]') + ' r=%f, p-value=%f' % pr
		print blue('[Spearman correlation]') + ' r=%f, p-value=%f' % sr
	# label = 
	plt.figure()
	plt.plot(v1, v2, 'o')
	plt.ylabel(var2)
	plt.xlabel(var1)
	plt.tight_layout()
	plt.show()


def do_remove_rotation(system, prot=None, nrem=1, fwhm=False, rhk=False, fix_p=True, full=False):
	try:
		system.per
	except AttributeError:
		if fwhm:
			system.per = gls(system, quantity='fwhm')
		if rhk:
			system.per = gls(system, quantity='rhk')
		else:
			system.per = gls(system)

	if prot is None:
		prot = system.per.get_peaks(output_period=True)[1]
	else:
		prot = float(prot)

	quantity = 'RV'
	if fwhm: 
		system.vrad = system.extras.fwhm
		system.error = 2.35 * system.error
		quantity = 'FWHM'
	elif rhk:
		system.vrad = system.extras.rhk
		system.error = system.extras.sig_rhk
		quantity = 'RHK'

	if fix_p:
		msg = blue('INFO: ') + 'Removing Prot=%.2f days (exactly) from %s \n' % (prot, quantity)
	else:
		msg = blue('INFO: ') + 'Removing Prot=%.2f days (fitting it) from %s\n' % (prot, quantity)
	msg += blue('    : ') + 'plus %d harmonics' % (nrem-1)
	clogger.info(msg)

	vrad_mean = system.vrad.mean()
	t, v, err = system.time, system.vrad - vrad_mean, system.error # temporaries

	if nrem == 1:
		def func(par, *args):
			if fix_p:
				As, Ac = par
				P = prot
			else:
				As, Ac, P = par
			# P = args
			return (v - (As*np.sin(2.*pi*t/P) + Ac*np.cos(2.*pi*t/P)) ) / err

		starting_values = [0.001, 0.001]
		if not fix_p: starting_values = starting_values + [prot]
		rot_param = leastsq(func, starting_values, maxfev=50000)[0]
		print rot_param

	elif nrem == 2:	
		def func2(par, *args):
			if fix_p:
				As1, Ac1, As2, Ac2 = par
				P = prot
			else:
				As1, Ac1, As2, Ac2, P = par
			# P = args[0]
			Po2 = P/2.
			return (v - (As1*np.sin(2.*pi*t/P) + Ac1*np.cos(2.*pi*t/P)) - (As2*np.sin(2.*pi*t/Po2) + Ac2*np.cos(2.*pi*t/Po2)) ) / err

		starting_values = [0.1, 0.1, 0.01, 0.01]
		if not fix_p: starting_values = starting_values + [prot]
		rot_param = leastsq(func2, starting_values, maxfev=500000)[0]
		print rot_param

	if fix_p: 
		prot_fit = prot
	else: 
		prot_fit = rot_param[4] # the period that we ended up removing

	print 'Phase = ', np.arctan2(rot_param[1], rot_param[0])

	system.per._plot(verts=[prot_fit/i for i in range(1,nrem+1)], doFAP=True)

	plt.figure()
	plt.subplot(211)
	plt.plot(system.time, system.vrad, 'o')
	plt.subplot(212)
	xx = np.linspace(t.min(), t.max(), 500)

	if nrem==1:
		if fix_p:
			As1, Ac1 = rot_param
			P = prot
		else:
			As1, Ac1, P = rot_param
		
		plt.plot(xx, (As1*np.sin(2.*pi*xx/P) + Ac1*np.cos(2.*pi*xx/P)), 'o-')

	if nrem==2:
		if fix_p:
			As1, Ac1, As2, Ac2 = rot_param
			P = prot
		else:
			As1, Ac1, As2, Ac2, P = rot_param

		plt.plot(xx, (As1*np.sin(2.*pi*xx/P) + Ac1*np.cos(2.*pi*xx/P)) + \
		             (As2*np.sin(4.*pi*xx/P) + Ac2*np.cos(4.*pi*xx/P)), 'o-')
	
	plt.errorbar(system.time, system.vrad - vrad_mean, yerr=err, fmt='ro')

	# plt.figure()
	# plt.subplot(111)
	# plt.plot(system.time, system.vrad, 'o')
	if nrem == 1: vrad_new = func(rot_param, prot)*err + vrad_mean
	elif nrem == 2: vrad_new = func2(rot_param, prot)*err + vrad_mean
	# plt.plot(system.time, vrad_new, 'go')

	# plt.show()

	newsystem = copy.copy(system)
	newsystem.vrad = vrad_new
	per = gls(newsystem)
	per._plot(doFAP=True, verts=[prot/i for i in range(1,nrem+1)]+[18])

	msg = blue('INFO: ') + 'Setting new RV array'
	clogger.info(msg)
	system.vrad = newsystem.vrad

	# msg = blue('INFO: ') + 'Removing harmonics also from FWHM'
	# clogger.info(msg)

	# d = system.extras._asdict() # asdict because namedtuple is immutable
	# field = 'fwhm'
	# d[field] = system.extras.fwhm - system.extras.fwhm.mean() \
	#            - (As1*np.sin(2.*pi*t/P) + Ac1*np.cos(2.*pi*t/P)) + (As2*np.sin(4.*pi*t/P) + Ac2*np.cos(4.*pi*t/P))

	# extra = namedtuple('Extra', system.extras_names, verbose=False)
	# system.extras = extra(**d)


	if full:
		##### do a "global" fit with 1 keplerian plus removing rotation harmonics
		xx = np.linspace(t.min(), t.max(), 500)
		vel1 = zeros_like(xx) 
		vel = zeros_like(system.time)

		def func3(params, return_model=False, times=None):
			""" Fitness function for 1 planet model plus rotational harmonics """
			#print params
			P, K, ecc, omega, T0, gam, As1, Ac1, As2, Ac2, Prot = params
			Po2 = Prot / 2.
			# Po3 = Prot / 3.
			#print ecc
			if (ecc>1 or ecc<0): return 1e99
			# if any(e>1 or e<0 for e in ecc): return 1e99
			if return_model:
				get_rvn(times, P, K, ecc, omega, T0, gam, vel1)
				return vel1 + \
				        (As1*np.sin(2.*pi*times/Prot) + Ac1*np.cos(2.*pi*times/Prot)) + \
				        (As2*np.sin(2.*pi*times/Po2) + Ac2*np.cos(2.*pi*times/Po2)) #+ \
				        # (As3*np.sin(2.*pi*times/Po3) + Ac3*np.cos(2.*pi*times/Po3))
			else:
				get_rvn(t, P, K, ecc, omega, T0, gam, vel)
				return v - vel - \
				        (As1*np.sin(2.*pi*t/Prot) + Ac1*np.cos(2.*pi*t/Prot)) - \
				        (As2*np.sin(2.*pi*t/Po2) + Ac2*np.cos(2.*pi*t/Po2)) #- \
				        # (As3*np.sin(2.*pi*t/Po3) + Ac3*np.cos(2.*pi*t/Po3))

		# initial parameters [P, K, ecc, omega, T0, gam, As1, Ac1, As2, Ac2, As3, Ac3, Prot]
		x0 = [18., 0.005, 0.2, 0., t.min(), 0., 0.01, 0.01, 0.005, 0.005, prot]
		best_param = leastsq(func3, x0)[0]
		print best_param

		plt.figure()
		plt.subplot(211)
		# data
		plt.plot(system.time, system.vrad, 'o')
		# fit
		plt.plot(xx, func3(best_param, return_model=True, times=xx), 'g-')
		plt.ylim([-0.01, 0.015])
		# keplerian curve only
		plt.subplot(212)
		P, K, ecc, omega, T0, gam, As1, Ac1, As2, Ac2, Prot = best_param
		get_rvn(xx, P, K, ecc, omega, T0, gam, vel1)
		plt.plot(xx, vel1, 'k-', label='keplerian, P='+str(P))
		plt.legend()
		plt.ylim([-0.01, 0.015])
		plt.show()

		newsystem = copy.copy(system)
		newsystem.vrad = func3(best_param)
		per = gls(newsystem)
		per._plot()
		system.vrad = newsystem.vrad



def do_Dawson_Fabrycky(system):
	
	if not periodogram_DF_available:
		msg = red('ERROR: ') + 'This extension is not available. Something went wrong on install...'
		clogger.fatal(msg)
		return

	time, rv, err = system.time, system.vrad, system.error
	ofac = 2.0

	def specwindow(freq,time):
		""" Calculate the spectral window function and its phase angles """
		n = len(time)
		W = [sum([np.exp(-2.j*pi*f*t) for t in time])/float(n) for f in freq]
		amp = [np.sqrt(t.real*t.real + t.imag*t.imag) for t in W]
		phase = [np.arctan2(t.imag,t.real) for t in W]
		return amp,phase

	plow = 0.5
	n = len(err)

	### GET THE REAL PERIODOGRAM
	freq,power,a_cos,b_sin,c_cte,phi,fNy,xdif = periodogram_DF(time, rv, err, ofac, plow)
	freq_real, power_real = freq, power 
	### GET THE WINDOW FUNCTION AT THOSE FREQUENCIES + plot dials
	amp,phase = specwindow(freq,time)

	figure(num=1,figsize=(12,8))

	#### GET 3 Maximum peaks + create fake data
	nf = len(freq)
	temp = sorted(zip(power,freq,a_cos,b_sin,c_cte,phi))
	fmax1 = (temp[nf-1])[1]

	fmax2 = fmax1
	i = 1
	while abs(fmax2 - fmax1) < 0.01:
		i += 1
		fmax2 = (temp[nf-i])[1]

	fmax3 = fmax2
	j = i
	while abs(fmax3 - fmax2) < 0.01 or abs(fmax3 - fmax1) < 0.01:
		j += 1
		fmax3 = (temp[nf-j])[1]

	print 'Maximum frequencies at: ',fmax1,fmax2,fmax3
	print 'That is periods at: ',round(1./fmax1,2),round(1./fmax2,2),round(1./fmax3,2)

	a1 = (temp[nf-1])[2]
	a2 = (temp[nf-i])[2]
	a3 = (temp[nf-j])[2]
	b1 = (temp[nf-1])[3]
	b2 = (temp[nf-i])[3]
	b3 = (temp[nf-j])[3]
	c1 = (temp[nf-1])[4]
	c2 = (temp[nf-i])[4]
	c3 = (temp[nf-j])[4]
	ph1 = (temp[nf-1])[5]
	ph2 = (temp[nf-i])[5]
	ph3 = (temp[nf-j])[5]

	timemin=int(min(time))
	timemax=int(max(time))
	timefake=frange(timemin-10,timemax+10,0.05)
	timefake = time

	xdiff = max(time) - min(time)

	rv_fake1 = array([a1*cos(fmax1*2.*pi*i) + b1*sin(fmax1*2.*pi*i) + c1 for i in timefake])
	rv_fake2 = array([a2*cos(fmax2*2.*pi*i) + b2*sin(fmax2*2.*pi*i) + c2 for i in timefake])
	rv_fake3 = array([a3*cos(fmax3*2.*pi*i) + b3*sin(fmax3*2.*pi*i) + c3 for i in timefake])
	#errfake = [0.001 for i in timefake]
	errfake = err

	### PLOT REAL PERIODOGRAM + DIALS
	figure(num = 1)

	subplot(4,1,1)
	title('window function + periodogram')
	semilogx(1/freq,amp,'r-', alpha=0.3)
	semilogx(1/freq,power,'k-')

	semilogx(1./fmax1,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax2,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax3,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)

	semilogx([1./fmax1,1./fmax1+0.025*cos(ph1)],[max(power)+0.1,max(power)+0.1+0.025*sin(ph1)],'k-',lw=1)
	semilogx([1./fmax2,1./fmax2+0.025*cos(ph2)],[max(power)+0.1,max(power)+0.1+0.025*sin(ph2)],'k-',lw=1)
	semilogx([1./fmax3,1./fmax3+0.025*cos(ph3)],[max(power)+0.1,max(power)+0.1+0.025*sin(ph3)],'k-',lw=1)

	xlim(plow,xdiff*ofac)
	ylim(0.0,max(power)+0.2)

	### PLOT FAKE PERIODOGRAMS + DIALS
	#### 1st FAKE
	freq,power,a_cos,b_sin,c_cte,phi,fNy,xdif = periodogram_DF(timefake, rv_fake1, errfake, ofac, plow)

	ind1 = list(freq).index(fmax1)
	ind2 = list(freq).index(fmax2)
	ind3 = list(freq).index(fmax3)

	subplot(4,1,2)

	semilogx([1./i for i in freq],power,'k-')
	fill_between(1/freq_real, power_real, 0., color='k', alpha=0.5)

	semilogx(1./fmax1,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax2,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax3,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)

	semilogx([1./fmax1,1./fmax1+0.045*cos(phi[ind1])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind1])],'k-',lw=1)
	semilogx([1./fmax2,1./fmax2+0.045*cos(phi[ind2])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind2])],'k-',lw=1)
	semilogx([1./fmax3,1./fmax3+0.045*cos(phi[ind3])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind3])],'k-',lw=1)

	xlim(plow,xdiff*ofac)
	ylim(0.0,max(power)+0.2)

	#### 2nd FAKE
	freq,power,a_cos,b_sin,c_cte,phi,fNy,xdif = periodogram_DF(timefake, rv_fake2, errfake, ofac, plow)

	ind1 = list(freq).index(fmax1)
	ind2 = list(freq).index(fmax2)
	ind3 = list(freq).index(fmax3)

	subplot(4,1,3)

	semilogx([1./i for i in freq],power,'k-')

	semilogx(1./fmax1,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax2,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax3,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)

	semilogx([1./fmax1,1./fmax1+0.045*cos(phi[ind1])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind1])],'k-',lw=1)
	semilogx([1./fmax2,1./fmax2+0.045*cos(phi[ind2])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind2])],'k-',lw=1)
	semilogx([1./fmax3,1./fmax3+0.045*cos(phi[ind3])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind3])],'k-',lw=1)

	xlim(plow,xdiff*ofac)
	ylim(0.0,max(power)+0.2)

	#### 3rd FAKE
	freq,power,a_cos,b_sin,c_cte,phi,fNy,xdif = periodogram_DF(timefake, rv_fake3, errfake, ofac, plow)

	ind1 = list(freq).index(fmax1)
	ind2 = list(freq).index(fmax2)
	ind3 = list(freq).index(fmax3)

	subplot(4,1,4)

	semilogx([1./i for i in freq],power,'k-')

	semilogx(1./fmax1,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax2,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax3,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)

	semilogx([1./fmax1,1./fmax1+0.045*cos(phi[ind1])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind1])],'k-',lw=1)
	semilogx([1./fmax2,1./fmax2+0.045*cos(phi[ind2])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind2])],'k-',lw=1)
	semilogx([1./fmax3,1./fmax3+0.045*cos(phi[ind3])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind3])],'k-',lw=1)

	xlim(plow,xdiff*ofac)
	ylim(0.0,max(power)+0.2)

	# savefig(name+'_DF.ps',orientation = 'Landscape')
	show()


def do_create_planets():
	msg = blue('INFO: ') + 'Starting the planet creator\n'
	clogger.info(msg)

	try:
		nplanets = int(raw_input('How many planets? '))
	except ValueError:
		msg = red('ERROR: ') + "I don't know how many planets to make\n"
		clogger.fatal(msg)
		return
	

	if ask_yes_no('Specify the periods (y/N)? ', default=False):
	    periods = []
	    for i in range(nplanets):
	        p = float(raw_input('\tperiod %d:  '%(i+1)))
	        periods.append(p)
	else: 
	    periods = None
	
	if ask_yes_no('Specify the eccentricities (y/N)? ', default=False):
	    eccentricities = []
	    for i in range(nplanets):
	        e = float(raw_input('\teccentricity %d:  '%(i+1)))
	        eccentricities.append(e)
	else: 
	    eccentricities = None
	        
	if ask_yes_no('Use sampling from file (y/N)? ', default=False):
		filename = raw_input('Which file? ')
		with open(filename) as f:
			nobs = len(f.readlines())
	else:
		filename = None
		obs = raw_input('How many observed RVs? ("r" for random) ')
		nobs = np.random.randint(30, 220) if obs in ('r','') else int(obs)

	noise = 'wrong_option'
	while noise not in ('w', '', 'r', 'R', 'wr', 'rw'):
		noise = raw_input('What type of noise? (W / r / wr) ')
		if noise in ('w', ''): 
		    type_noise = 'white'
		elif noise in ('r', 'R'): 
		    type_noise = 'correlated'
		elif noise in ('wr', 'rw'): 
		    type_noise = 'white + correlated'
		else: 
		    print "I don't understand that. Try 'w', 'r' or 'wr'."

	print 
	msg =  blue('INFO: ') + 'Generating %d planet(s)\n' % (nplanets)
	if filename is not None:
		msg += blue('    : ') + '-> sampling from %s\n' % (filename)	
	else:
		msg += blue('    : ') + '-> randomly spaced sampling\n'
	msg += blue('    : ') + '-> %d observations\n' % (nobs)
	msg += blue('    : ') + '-> %s noise\n' % (type_noise)
	clogger.info(msg)	

	# get rv curve with random parameters
	def gen_rv(nplanets, nobs, sampling_file, periods, eccentricities, type_noise, temp):
	    
	    if sampling_file is not None: # the sampling is set in the file
	        times_sampled_from_file = np.loadtxt(sampling_file, usecols=(0,), skiprows=1)
	        times = times_sampled_from_file
	        times_full = np.linspace(min(times_sampled_from_file), max(times_sampled_from_file), 1000)
	    else: # the sampling is random
	        times_full = np.linspace(2449460, 2452860, 1000)
	        # sample (N points) randomly from the data to produce unevenly sampled time series
	        rand_args = [ i for i in sorted(random.sample(xrange(len(times_full)), nobs)) ]
	        times = times_full[rand_args]

	    vel_full = zeros_like(times_full) # rv at oversampled times
	    vel_full_each = zeros_like(times_full) # rv at oversampled times for each planet
	    vel_all = np.zeros((len(times), nplanets)) # rv at sample times for each planet
	    vel_each = zeros_like(times) # rv at sample times for each planet (temporary)
	    vel_total = zeros_like(times) # rv at sample times

	    output = '# planet %d with P=%f, K=%f, e=%f, w=%f, t0=%f\n'
	    for planet in range(nplanets):
	        if periods is None: 
	            P = np.random.rand()*998. + 2. # random period, U(2, 998)
	        else: 
	            P = periods[planet] # user-provided period
	        
	        if eccentricities is None:
	        	ecc = np.random.beta(0.867, 3.03)  # from Kipping 2013
	        else:
	        	ecc = eccentricities[planet] # user-provided eccentricity

	        K = np.random.rand()*50. + 2.

	        omega = np.random.rand()*2.*pi 

	        if type_noise is 'white':
	            noise = np.random.randn(len(vel_each))*3. # sigma??????

	        # log
	        temp.write(output % (planet+1, P, K, ecc, omega, 2452000))

	        # get RVs for this planet
	        get_rvn(times, P, K, ecc, omega, 2452000, 0., vel_each)
	        vel_each = vel_each + noise # add noise

	        # get RVs for this planet at oversampled times
	        get_rvn(times_full, P, K, ecc, omega, 2452000, 0., vel_full_each)

	        # store
	        vel_all[:,planet] = vel_each
	        vel_total += vel_each
	        vel_full += vel_full_each

	    
	    return times_full, vel_full, times, vel_total, vel_all, noise

	# wrapper function that takes care of adding noise and logging to file
	def generate_planets(nplanets, nobs, sampling_file, type_noise, periods, eccentricities):
	    # create temporary file to store simulation info
	    tf = tempfile.NamedTemporaryFile(dir='./')
	    tf.write('## Planets created on ')
	    tf.write(time.strftime("%d/%m/%Y - %H:%M:%S\n#\n"))
	    
	    times_full, vel_full, times, vel, vel_all, noise = \
	                gen_rv(nplanets, nobs, sampling_file, periods, eccentricities, type_noise, tf)
	    
	    tf.write('#\n## %d observations, %s noise\n#\n'%(nobs, type_noise))
	    tf.write('# bjd \t vrad(km/s) \t svrad(km/s)\n')

	    for t, v, e in zip(times, vel, noise):
	        tf.write('%f\t%f\t%f\n' % (t, v*1e-3, abs(e*1e-3)))
	    tf.delete = False
	    print 'Output to %s' % tf.name
	    tf.close()
	        
	    return times_full, vel_full, times, vel, vel_all


	times_full, vel_full, times, vel1, vel2 = generate_planets(nplanets, nobs, filename, type_noise, periods, eccentricities)


def load_plugin(plugin):
	from runpy import run_module
	current_dir = os.path.dirname(__file__)

	if plugin in os.listdir(current_dir+'/plugins'):
		msg = blue('INFO: ') + 'Loading plugin "%s"\n' % plugin
		clogger.info(msg)

		relative_module = 'OPEN.plugins.' + plugin + '.run'
		run_module(relative_module)
		# module = __import__(relative_module, fromlist=['OPEN.plugins'])
		# func = getattr(module, plugin)
		# func()