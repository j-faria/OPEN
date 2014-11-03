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

# see http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
from scipy.optimize import leastsq
from scipy.stats.stats import spearmanr, pearsonr
from scipy.stats import nanmean, nanstd
from scipy.odr import odrpack

# intra-package imports
from .classes import MCMC_dream, MCMC_nest
from logger import clogger
from ext.get_rvN import get_rvn
from .periodograms import ls_PressRybicki, gls
try:
	from ext.periodogram_DF import periodogram_DF
	periodogram_DF_available = True
except ImportError:
	periodogram_DF_available = False
try:
	from ext.periodogram_CLEAN import clean
	periodogram_CLEAN_available = True
except ImportError:
	periodogram_CLEAN_available = False

from shell_colors import yellow, red, blue
from .utils import julian_day_to_date, ask_yes_no, get_number_cores

timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
pi = np.pi

def do_fit(system, verbose):
	try:
		degree = system.model['d']
		## trend removal
		# did we already remove this trend?
		if (system.model.has_key('drift') and degree == len(system.model['drift'])-1):
			msg = yellow('RESULT: ') + 'Fit of degree %d done! Coefficients:' % degree
			clogger.info(msg)
			msg = yellow('      : ') + str(system.model['drift'])
			clogger.info(msg)
			return

		# did we already remove a different trend?
		if (system.model.has_key('drift') and degree != len(system.model['drift'])-1):
			system.vrad = system.vrad_full

		if degree > 0:
			z = np.polyfit(system.time, system.vrad, degree)
			# save info to system
			system.model['drift'] = z
			# plot fitted drift
			system.do_plot_drift()
			# change system RVs accordingly
			poly = np.poly1d(z)
			system.vrad = system.vrad - poly(system.time)
			## output some information
			msg = yellow('RESULT: ') + 'Fit of degree %d done! Coefficients:' % degree
			clogger.info(msg)
			msg = yellow('      : ') + str(z)
			clogger.info(msg)
	except TypeError:
		msg = yellow('Warning: ') + 'To remove trends, run mod before fit.\n'
		clogger.error(msg)

	kep = system.model['k']
	if (kep == 0): 
		return

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
		# system.per = ls_PressRybicki(system)  # this periodogram is very fast
		system.per = gls(system)  # this periodogram is very fast
	finally:
		peak = system.per.get_peaks(output_period=True)[0]

	peak_val = system.per.power.max()
	peak_fap = system.per.FAP(peak_val)
	if peak_fap < 0.01:
		msg = blue('INFO: ') + 'Peak found at %f days with FAP ~ %e' % (peak[0], peak_fap)
		clogger.info(msg)
	else:
		msg = yellow('Warning: ') + 'No peaks found with FAP < 1%'
		clogger.error(msg)
		return None

	p0 = peak[0] # initial guess period
	k0 = system.vrad.max()
	e0 = 0.0
	om0 = 0.0
	t0 = system.time.min()
	g0 = system.vrad.mean()

	msg = blue('INFO: ') + 'Setting initial guesses...'  #, p0, k0, e0, om0, t0, g0
	clogger.info(msg)
	initial_guess = [p0, k0, e0, om0, t0, g0]

	## call levenberg markardt fit
	lm = do_lm(system, initial_guess)
	lm_par = lm[0]
	# lm_par[1] = abs(lm_par[1])
	# lm_par[3] = np.rad2deg(lm_par[3])

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
	par_all = lm_par  # backup

	# j=2
	# while j<3:
	# 	# is there a relevant peak in the periodogram of the residuals?
	# 	system.per_resid = gls(system, quantity='resid')
	# 	# system.per_resid._plot()
	# 	peak = system.per_resid.get_peaks(output_period=True)[0]
	# 	peak_val = system.per_resid.power.max()
	# 	peak_fap = system.per_resid.FAP(peak_val)
	# 	if peak_fap < 0.01:
	# 		msg = blue('INFO: ') + 'Peak found at %f days with FAP ~ %e' % (peak[0], peak_fap)
	# 		clogger.info(msg)
	# 	else:
	# 		msg = yellow('Warning: ') + 'No peaks found in the periodogram of the residuals with FAP < 1%'
	# 		clogger.error(msg)
	# 		break

	# 	# 	new_system = copy.deepcopy(system)
	# 	# 	new_system.vrad = new_system.fit['residuals']

	# 	p0 = 5000 #peak[0] # initial guess period
	# 	k0 = 370.
	# 	e0 = 0.8
	# 	om0 = 300.
	# 	t0 = system.time.mean()
	# 	g0 = system.vrad.mean()
	# 	msg = blue('INFO: ') + 'Setting initial guesses for another planet...' #, p0, k0, e0, om0, t0, g0
	# 	initial_guess = np.append(par_all, [p0, k0, e0, om0, t0, g0])
		
	# 	print initial_guess
	# 	## call levenberg markardt fit
	# 	lm = do_lm(system, initial_guess)
	# 	lm_par = lm[0]
		
	# 	## loop over planets
	# 	msg = yellow('RESULT: ') + 'Best fit is'
	# 	clogger.info(msg)
	# 	print("%3s %12s %10s %10s %10s %15s %9s" % \
	# 		('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
	# 	for i, planet in enumerate(list(ascii_lowercase)[:j]):
	# 		P, K, ecc, omega, T0, gam = [lm_par[k::6] for k in range(6)]
	# 		print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )

	# 	par_all = np.append(par_all, lm_par)
	# 	# save fit in the system
	# 	system.save_fit(par_all, 0.)
	# 	j+=1


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

def do_demc(system, zfile, burnin=0):
	import demc.dream as dream
	from time import sleep
	from re import findall

	if zfile is not None:
		msg = blue('INFO: ') + 'Loading chains from file ' + zfile
		clogger.info(msg)
		results = MCMC_dream(zfile, burnin=burnin)
		results.print_best_solution()
		results.do_plot_some_solutions(system)
		return results

	try:
		degree = system.model['d']
		keplerians = system.model['k']
	except TypeError:
		msg = red('Error: ') + 'Need to run mod before demc. '
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

	results = MCMC_dream('OPEN/demc/problem1_chain*', burnin=burnin)

	## output best solution
	results.print_best_solution()

	print  # newline
	msg = blue('INFO: ') + 'Estimating kernels...'
	clogger.info(msg)

	# p = results.kde1d('period_0', npoints=300)
	# ind = np.linspace(0.2, 5000, 300)

	results.do_plot_best(system)

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
	from deap import base, creator, tools, algorithms
	
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

	def chi2_n_leastsq(params):
		""" Fitness function for N planet model """
		P, K, ecc, omega, T0, gam = [params[i::6] for i in range(6)]

		if any(e > 1 or e < 0 for e in ecc): return 1e99
		# if any(k < 0 for k in K): return 1e99

		get_rvn(system.time, P, K, ecc, omega, T0, gam[0], vel)
		return (system.vrad - vel) / system.error

	def chi2_2(params):
		P, K, ecc, omega, T0, gam = params
		get_rvn(system.time, P, K, ecc, omega, T0, gam, vel)
		return system.vrad - vel

	return leastsq(chi2_n_leastsq, x0, full_output=0, ftol=1e-15, maxfev=int(1e6))


def do_multinest(system, user, gp, resume=False, ncpu=None, training=None, lin=None, doplot=True, feed=False):
	from time import sleep, time

	def get_multinest_output(root, nplanets, context='111'):
		msg = blue('INFO: ') + 'Analysing output...'
		clogger.info(msg)

		with open(root+'stats.dat') as f:
			stats = f.readlines()
			nlines = len(stats)

		npar = (nlines - 10) / 3
		if context[2] == '3': return

		try:
			NS_lnE = float(stats[0].split()[-3])
			NS_lnE_error = float(stats[0].split()[-1])

			INS_lnE = float(stats[1].split()[-3])
			INS_lnE_error = float(stats[1].split()[-1])
		except ValueError:
			NS_lnE, NS_lnE_error = 0, 0
			INS_lnE, INS_lnE_error = 0, 0

		## mean params
		par_mean = [float(s.split()[1]) for s in stats[4:4+npar]]
		# P_mean, K_mean, ecc_mean, omega_mean, t0_mean, vsys_mean = par_mean
		par_sigma = [float(s.split()[2]) for s in stats[4:4+npar]]
		# P_sigma, K_sigma, ecc_sigma, omega_sigma, t0_sigma, vsys_sigma = par_sigma

		## MLE
		start, end = 4+npar+3, 4+npar+3+npar
		par_mle = [float(s.split()[1]) for s in stats[start:end]]
		# P_mle, K_mle, ecc_mle, omega_mle, t0_mle, vsys_mle = par_mle

		## MAP
		start, end = 4+2*3+2*npar, 4+2*3+3*npar
		par_map = [float(s.split()[1]) for s in stats[start:end]]
		# P_map, K_map, ecc_map, omega_map, t0_map, vsys_map = par_map

		msg = yellow('RESULT: ') + 'Parameters summary'
		clogger.info(msg)

		print '%8s %14s %9s %14s %14s' % ('', 'mean', '+- sigma', 'ML', 'MAP')
		## loop over planets
		i = -1
		for i, planet in enumerate(list(ascii_lowercase)[:nplanets]):
			print yellow(planet)
			print '%8s %14.3f %9.3f %14.3f %14.3f' % ('P',     par_mean[5*i], par_sigma[5*i], par_mle[5*i], par_map[5*i])
			print '%8s %14.3f %9.3f %14.3f %14.3f' % ('K',     par_mean[5*i+1], par_sigma[5*i+1], par_mle[5*i+1], par_map[5*i+1])
			print '%8s %14.3f %9.3f %14.3f %14.3f' % ('ecc',   par_mean[5*i+2], par_sigma[5*i+2], par_mle[5*i+2], par_map[5*i+2])
			print '%8s %14.3f %9.3f %14.3f %14.3f' % ('omega', par_mean[5*i+3], par_sigma[5*i+3], par_mle[5*i+3], par_map[5*i+3])
			print '%8s %14.2f %9.2f %14.2f %14.2f' % ('t0',    par_mean[5*i+4], par_sigma[5*i+4], par_mle[5*i+4], par_map[5*i+4])
		print yellow('system')
		if context[2] == '2':
			# jitter parameter
			print '%8s %14.3f %9.3f %14.3f %14.3f' % ('jitter', par_mean[-2], par_sigma[-2], par_mle[-2], par_map[-2])
		if context[0] == '1':
			# in this case, the vsys parameters is the last one
			print '%8s %14.3f %9.3f %14.3f %14.3f' % ('vsys', par_mean[-1], par_sigma[-1], par_mle[-1], par_map[-1])
		elif context[0] == '2':
			# in this case, the vsys is before the hyperparameters
			print '%8s %14.3f %9.3f %14.3f %14.3f' % ('vsys', par_mean[5*i+5], par_sigma[5*i+5], par_mle[5*i+5], par_map[5*i+5])

	available_cpu = get_number_cores()
	if (ncpu is None) or (ncpu > available_cpu): 
		ncpu = available_cpu

	msg = blue('INFO: ') + 'Transfering data to MultiNest...'
	clogger.info(msg)

	# training holds the quantity on which to train the GP beforehand or None
	if training:
		msg = blue('    : ') + 'Will train Gaussian Process on %s...' % training
		clogger.info(msg)
	# including linear dependence in the model
	if lin:
		msg = blue('    : ') + 'Model includes linear dependence on %s...' % lin
		clogger.info(msg)

	# write data to file to be read by MultiNest
	nest_filename = 'input.rv'
	nest_header = 'file automatically generated for MultiNest analysis, ' + timestamp

	d = system.provenance
	nest_header += '\n' + str(len(d))  # number of files (observatories)
	# this is a hack otherwise the dict values are not ordered the right way
	sizes_each_file = [d[k][0] if (d[k][1]==0) else (d[k][0]-d[k][1]) for k in sorted(d.keys())]
	sizes = len(sizes_each_file) * '%d ' % tuple(sizes_each_file)
	nest_header += '\n' + sizes  # number measurements in each file
	nest_header += '\n' + str(len(system.time))  # total number measurements

	if training:
		i = system.extras._fields.index(training) # index corresponding to the quantity
		savetxt(nest_filename, zip(system.time, system.vrad, system.error, system.extras[i]),
				header=nest_header,
				fmt=['%12.6f', '%7.5f', '%7.5f', '%7.5f'])
	elif lin:
		lvars = lin.split(',')
		i = [system.extras._fields.index(l) for l in lvars]
		if len(i) == 1:
			v1 = [system.extras[ii] for ii in i][0]
			out = zip(system.time, system.vrad, system.error, v1)
		elif len(i) == 2:
			v1, v2 = [system.extras[ii] for ii in i]
			out = zip(system.time, system.vrad, system.error, v1, v2)
		elif len(i) == 3:
			v1, v2, v3 = [system.extras[ii] for ii in i]
			out = zip(system.time, system.vrad, system.error, v1, v2, v3)

		savetxt(nest_filename, out,
				header=nest_header,
				fmt=['%12.6f', '%7.5f', '%7.5f'] + ['%7.5f']*len(i))
	else:
		savetxt(nest_filename, zip(system.time, system.vrad, system.error),
				header=nest_header,
				fmt=['%12.6f', '%7.5f', '%7.5f'])

	try:
		del system.results
	except AttributeError:
		pass

	if user:
		# user is controlling and editing the namelist, we just start 
		# multinest once with whatever is in there and read the output

		# but first we need some info from the namelist
		with open('OPEN/multinest/namelist1') as f:
			namelist_lines = f.readlines()

		# the root of output files which may have changed
		l1 = [line for line in namelist_lines if 'nest_root' in line][0]
		# split the assignment, strip of newline, strip of ', strip of "
		root_path = l1.split('=')[1].strip().strip("'").strip('"')

		# whether or not the model is a GP
		l1 = [line for line in namelist_lines 
		             if ('nest_context' in line) and not line.strip().startswith('!')][0]
		user_gp_context = int(l1.split('=')[1].strip()[0])

		# how many planets are in the model
		nplanets = int(l1.split('=')[1].strip()[1])

		# the full context
		full_context = l1.split('=')[1].strip()

		# if resuming from a previous run, set the appropriate namelist flag
		if resume:
			replacer = '    nest_resume=.true.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_resume' in line: print replacer,
				else: print line,
		else:
			replacer = '    nest_resume=.false.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_resume' in line: print replacer,
				else: print line,

		# if training the GP before, set the appropriate namelist flags
		if training:
			replacer1 = '    training = .true.\n'
			replacer2 = '    train_variable = \'%s\'\n' % training
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'training' in line: print replacer1,
				elif 'train_variable' in line: print replacer2,
				else: print line,
		else:
			replacer1 = '    training = .false.\n'
			replacer2 = '    train_variable = \'NONE\'\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'training' in line: print replacer1,
				elif 'train_variable' in line: print replacer2,
				else: print line,



		msg = blue('    : ') + 'Starting MultiNest (%d threads) ...' % (ncpu,)
		clogger.info(msg)

		start = time()
		cmd = 'mpirun -np %d ./OPEN/multinest/nest' % (ncpu,)
		rc = subprocess.call(cmd, shell=True)

		if (rc == 1): 
			msg = red('ERROR: ') + 'MultiNest terminated prematurely'
			clogger.fatal(msg)
			return

		print  # newline
		msg = blue('INFO: ') + 'MultiNest took %f s' % (time()-start)
		clogger.info(msg)

		# return
		get_multinest_output(root_path, nplanets, context=full_context)
		system.results = MCMC_nest(root_path, context=full_context)

		if doplot:
			system.results.do_plot_map(system)
			#system.results.do_plot_map_phased(system)

		# save fit in the system
		system.results.save_fit_to(system)


	else:
		# OPEN is in control and will try to run a full model selection
		# analysis, editing the namelist accordingly.

		nplanets = 0
		if gp: 
			context = 201
		else: 
			context = 101
		root = 'chains/nest-noplanet-'
		## automatically fill the necessary parameters in the inlist
		replacer1 = '    nest_context = %d\n' % context
		replacer2 = '    nest_root=\'%s\'\n' % root
		for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
			if ('nest_context' in line) and not line.strip().startswith('!'): print replacer1,
			elif 'nest_root' in line: print replacer2,
			elif (resume and 'nest_resume' in line): print '    nest_resume=.true.\n',
			else: print line,

		# if resuming from a previous run, set the appropriate namelist flag
		if resume:
			replacer = '    nest_resume=.true.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_resume' in line: print replacer,
				else: print line,
		else:
			replacer = '    nest_resume=.false.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_resume' in line: print replacer,
				else: print line,

		# by default, feedback on the sampling progress is omitted, but
		# the user can still request it
		if feed:
			replacer = '    nest_fb=.true.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_fb' in line: print replacer,
				else: print line,
		else:
			replacer = '    nest_fb=.false.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_fb' in line: print replacer,
				else: print line,


		sleep(1)
		start = time()
		cmd = 'mpirun -np %d ./OPEN/multinest/nest' % (ncpu,)
		rc = subprocess.call(cmd, shell=True)

		if (rc == 1): 
			msg = red('ERROR: ') + 'MultiNest terminated prematurely'
			clogger.fatal(msg)
			return

		print  # newline
		msg = blue('INFO: ') + 'MultiNest took %f s' % (time()-start)
		clogger.info(msg)

		get_multinest_output(root, nplanets, context=str(context))
		results_constant = MCMC_nest(root, context=str(context))
		if doplot:
			results_constant.do_plot_map(system)

		constant_lnE = results_constant.NS_lnE

		##############################################################################
		print
		msg = blue('    : ') + 'Starting MultiNest for 1-planet model.\n'
		msg += blue('    : ') + '(using %d threads). Please wait...' % (ncpu,)
		clogger.info(msg)

		nplanets = 1
		if gp: 
			context = 211
		else: 
			context = 111
		root = 'chains/nest-1planet-'
		## automatically fill the necessary parameters in the inlist
		replacer1 = '    nest_context = %d\n' % context
		replacer2 = '    nest_root=\'%s\'\n' % root
		for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
			if ('nest_context' in line) and not line.strip().startswith('!'): print replacer1,
			elif 'nest_root' in line: print replacer2,
			elif (resume and 'nest_resume' in line): print '    nest_resume=.true.\n',
			else: print line,

		# if resuming from a previous run, set the appropriate namelist flag
		if resume:
			replacer = '    nest_resume=.true.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_resume' in line: print replacer,
				else: print line,
		else:
			replacer = '    nest_resume=.false.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_resume' in line: print replacer,
				else: print line,

		# by default, feedback on the sampling progress is omitted, but
		# the user can still request it
		if feed:
			replacer = '    nest_fb=.true.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_fb' in line: print replacer,
				else: print line,
		else:
			replacer = '    nest_fb=.false.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_fb' in line: print replacer,
				else: print line,


		sleep(1)
		start = time()
		cmd = 'mpirun -np %d ./OPEN/multinest/nest' % (ncpu,)
		rc = subprocess.call(cmd, shell=True)

		if (rc == 1): 
			msg = red('ERROR: ') + 'MultiNest terminated prematurely'
			clogger.fatal(msg)
			return

		print  # newline
		msg = blue('INFO: ') + 'MultiNest took %f s' % (time()-start)
		clogger.info(msg)

		get_multinest_output(root, nplanets, context=str(context))
		results_one_planet = MCMC_nest(root, context=str(context))
		if doplot:
			results_one_planet.do_plot_map(system)

		one_planet_lnE = results_one_planet.NS_lnE

		##############################################################################
		print  # newline
		msg = blue('INFO: ') + 'Starting MultiNest for 2-planet model.\n'
		msg += blue('    : ') + '(using %d threads). Please wait...' % (ncpu,)
		clogger.info(msg)

		nplanets = 2
		if gp: 
			context = 221
		else: 
			context = 121
		root = 'chains/nest-2planet-'
		## automatically fill the necessary parameters in the inlist
		replacer1 = '    nest_context = %d\n' % context
		replacer2 = '    nest_root=\'%s\'\n' % root
		for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
			if ('nest_context' in line) and not line.strip().startswith('!'): print replacer1,
			elif 'nest_root' in line: print replacer2,
			elif (resume and 'nest_resume' in line): print '    nest_resume=.true.\n',
			else: print line,

		# if resuming from a previous run, set the appropriate namelist flag
		if resume:
			replacer = '    nest_resume=.true.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_resume' in line: print replacer,
				else: print line,
		else:
			replacer = '    nest_resume=.false.\n'
			for line in fileinput.input('OPEN/multinest/namelist1', inplace=True):
				if 'nest_resume' in line: print replacer,
				else: print line,

		sleep(1)
		start = time()
		cmd = 'mpirun -np %d ./OPEN/multinest/nest' % (ncpu,)
		rc = subprocess.call(cmd, shell=True)

		if (rc == 1): 
			msg = red('ERROR: ') + 'MultiNest terminated prematurely'
			clogger.fatal(msg)
			return

		print  # newline
		msg = blue('INFO: ') + 'MultiNest took %f s' % (time()-start)
		clogger.info(msg)


		get_multinest_output(root, nplanets, context=str(context))
		results_two_planet = MCMC_nest(root, context=str(context))
		if doplot:
			results_two_planet.do_plot_map(system)

		two_planet_lnE = results_two_planet.NS_lnE

		##############################################################################
		print  # newline
		msg = yellow('RESULT: ') + 'Evidence results'
		clogger.info(msg)

		msg1 = yellow('      : ') + 'cte lnE'.rjust(12) + ' = %12.6f' % (constant_lnE)
		msg2 = yellow('      : ') + '1 planet lnE = %12.6f' % (one_planet_lnE)
		msg3 = yellow('      : ') + '2 planet lnE = %12.6f' % (two_planet_lnE)

		## odds ratio
		import warnings
		with warnings.catch_warnings():
			warnings.filterwarnings('error')
			try:
				O11 = 1
				O21 = np.exp(one_planet_lnE) / np.exp(constant_lnE)
				O31 = np.exp(two_planet_lnE) / np.exp(constant_lnE)
				K = np.exp(one_planet_lnE) / np.exp(two_planet_lnE)
			except RuntimeWarning:
				try:
					import mpmath
				except ImportError:
					try:
						from sympy import mpmath
					except ImportError:
						warn = red('Warning: ') + 'Cannot exponentiate lnE; mpmath is not available'
						clogger.warning(warn)
						O11 = O21 = O31 = np.nan
				else:
					O11 = 1
					O21 = mpmath.exp(one_planet_lnE) / mpmath.exp(constant_lnE)
					O31 = mpmath.exp(two_planet_lnE) / mpmath.exp(constant_lnE)
					K = mpmath.exp(one_planet_lnE) / mpmath.exp(two_planet_lnE)
			finally:
				msg1 += ' '*5 + 'p = %-13.8f\n' % (O11/(O11+O21+O31), )
				msg2 += ' '*5 + 'p = %-13.8f\n' % (O21/(O11+O21+O31), )
				msg3 += ' '*5 + 'p = %-13.8f\n' % (O31/(O11+O21+O31), )
				msg = msg1 + msg2 + msg3
				clogger.info(msg)

		

		# clogger.info(msg)

		# msg =  yellow('      : ') + 'Odds ratio = %f \n' % (K)
		# clogger.info(msg)

		# jeffreys_scale = {(0, 1) : 'Negative', 
		#                   (1, 3) : 'Barely worth mentioning', 
		#                   (3, 10) : 'Substantial', 
		#                   (10, 30) :'Strong', 
		#                   (30, 100) : 'Very strong', 
		#                   (100, 1e99) : 'Decisive'}
		# kass_raftery_scale = {(0, 1) : 'Negative', 
		#                       (1, 3) : 'Not worth more than a bare mention', 
		#                   	  (3, 20) : 'Positive', 
		#                   	  (20, 150) :'Strong', 
		#                   	  (150, 1e99) : 'Very strong'}

		# print  # newline
		# msg =  yellow('      : ') + '%-15s\t%s' % ('Scale', 'strength of evidence supporting 1 planet')
		# clogger.info(msg)

		# msg =  yellow('      : ')
		# for key in jeffreys_scale:
		# 	if key[0] < K <= key[1]: 
		# 		msg += '%-15s\t%s' % ('Jeffreys', jeffreys_scale[key])
		# clogger.info(msg)

		# msg =  yellow('      : ')
		# for key in kass_raftery_scale:
		# 	if key[0] < K <= key[1]: 
		# 		msg += '%-15s\t%s' % ('Kass & Raftery', kass_raftery_scale[key])
		# clogger.info(msg)

		## save the fit with highest evidence to the system
		fits = {constant_lnE: results_constant,
		        one_planet_lnE: results_one_planet,
		        two_planet_lnE: results_two_planet}
		system.results = fits[sorted(fits, reverse=True)[0]]


	return

def do_correlate(system, vars=(), verbose=False, remove=False):
	"""
	Correlations between radial velocities and/or other diagnostics
	"""
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

	e1 = np.zeros_like(system.vrad)
	e2 = np.zeros_like(system.vrad)

	if var1 == 'vrad': 
		v1 = system.vrad
		e1 = system.error
	else:
		i = system.extras._fields.index(var1)
		v1 = system.extras[i]
		if var1 == 'rhk': e1 = system.extras.sig_rhk

	if var2 == 'vrad': 
		v2 = system.vrad
		e2 = system.error
	else:
		i = system.extras._fields.index(var2)
		v2 = system.extras[i]
		if var2 == 'rhk': e2 = system.extras.sig_rhk

	pr = pearsonr(v1, v2)
	sr = spearmanr(v1, v2)

	if verbose:
		print blue('[Pearson correlation]') + ' r=%f, p-value=%f' % pr
		print blue('[Spearman correlation]') + ' r=%f, p-value=%f' % sr

	plt.figure()
	plt.errorbar(v1, v2, xerr=e1, yerr=e2, fmt='o')
	if verbose:
		# non-weighted fit, OLS
		m, b = np.polyfit(v1, v2, 1)
		yp = np.polyval([m, b], v1)
		plt.plot(v1, yp, '-k')

		# weghted fit, ODR (only do if errors on variables)
		if not (np.count_nonzero(e1) == 0 or np.count_nonzero(e2) == 0):
			def f(B, x): return B[0]*x + B[1]
			linear = odrpack.Model(f)
			data = odrpack.Data(v1, v2, wd=1./e1**2, we=1./e2**2)
			odr = odrpack.ODR(data, linear, beta0=[1., 1.])
			output = odr.run()
			yp = np.polyval(output.beta, v1)
			plt.plot(v1, yp, '--k')

	plt.ylabel(var2)
	plt.xlabel(var1)
	plt.tight_layout()
	plt.show()

	if remove:
		assert (var1 == 'vrad' or var2 == 'vrad')
		if var1 == 'vrad':
			# switch up the arrays to make sure we remove things the right way
			# v2 will always be the RVs
			v1, v2 = v2, v1
			e1, e2 = e2, e1
		vs = [var1, var2]
		# index of the variable which is not vrad
		i = int(not vs.index('vrad'))

		msg = blue('INFO: ') + 'Removing linear dependence RV ~ %s' % vs[i]
		clogger.info(msg)

		# weghted fit, ODR (only do if errors on variables)
		if not (np.count_nonzero(e1) == 0 or np.count_nonzero(e2) == 0):
			def f(B, x): return B[0]*x + B[1]
			linear = odrpack.Model(f)
			data = odrpack.Data(v1, v2, wd=1./e1**2, we=1./e2**2)
			odr = odrpack.ODR(data, linear, beta0=[1., 1.])
			output = odr.run()
			yp = np.polyval(output.beta, v1)

			system.vrad = system.vrad - yp

		else:  # non-weighted fit, OLS
			m, b = np.polyfit(v1, v2, 1)
			yp = np.polyval([m, b], v1)

			system.vrad = system.vrad - yp



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

	elif nrem == 3:
		def func3(par, *args):
			if fix_p:
				As1, Ac1, As2, Ac2, As3, Ac3 = par
				P = prot
			else:
				As1, Ac1, As2, Ac2, P = par
			# P = args[0]
			Po2 = P/2.
			Po3 = P/3.
			return (v - (As1*np.sin(2.*pi*t/P) + Ac1*np.cos(2.*pi*t/P)) -
				        (As2*np.sin(2.*pi*t/Po2) + Ac2*np.cos(2.*pi*t/Po2)) -
				        (As3*np.sin(2.*pi*t/Po3) + Ac3*np.cos(2.*pi*t/Po3)) ) / err

		starting_values = [0.1, 0.1, 0.01, 0.01, 0.01, 0.01]
		if not fix_p: starting_values = starting_values + [prot]
		rot_param = leastsq(func3, starting_values, maxfev=500000)[0]
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
	
	if nrem==3:
		if fix_p:
			As1, Ac1, As2, Ac2, As3, Ac3 = rot_param
			P = prot
		else:
			As1, Ac1, As2, Ac2, As3, Ac3, P = rot_param

		plt.plot(xx, (As1*np.sin(2.*pi*xx/P) + Ac1*np.cos(2.*pi*xx/P)) + \
					 (As2*np.sin(4.*pi*xx/P) + Ac2*np.cos(4.*pi*xx/P)) + \
					 (As3*np.sin(6.*pi*xx/P) + Ac3*np.cos(6.*pi*xx/P)), 'o-')

	plt.errorbar(system.time, system.vrad - vrad_mean, yerr=err, fmt='ro')

	# plt.figure()
	# plt.subplot(111)
	# plt.plot(system.time, system.vrad, 'o')
	if nrem == 1: vrad_new = func(rot_param, prot)*err + vrad_mean
	elif nrem == 2: vrad_new = func2(rot_param, prot)*err + vrad_mean
	elif nrem == 3: vrad_new = func3(rot_param, prot)*err + vrad_mean
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



def get_rotation_period(system):
	"""
	Calculate rotation period from the Noyes (1984) relation
	"""
	# try to find name of star automatically
	import re
	filename = system.provenance.keys()[0]
	regex = re.compile("HD[0-9]*")
	temp_star_name = regex.findall(filename)[0]
	if ask_yes_no('Is "%s" the star (Y/n)? ' % temp_star_name, default=True):
		pass
	else:
		temp_star_name = raw_input('Name of the star: ').upper()

	print  # newline
	msg = blue('INFO: ') + 'Calculating rotation period for %s' % temp_star_name
	clogger.info(msg)

	data_path = '/home/joao/phd/data/'
	filepath = data_path + 'metadata/sample_simbad.rdb'
	found_it = False
	with open(filepath) as f:
		for line in f:
			if temp_star_name in line:
				found_it = True
				star_info = line.strip().split()

	if not found_it:
		msg = red('ERROR: ') + 'Cannot find B-V information for %s' % temp_star_name
		clogger.fatal(msg)
		return

	msg = blue('    : ') + '%s, V=%s  B=%s' % tuple(star_info[1:])
	clogger.info(msg)

	Vmag = float(star_info[2])
	Bmag = float(star_info[3])

	x = 1 - (Bmag - Vmag)
	## Noyes (1984), Eq 4
	log_tau = (1.362 - 0.166*x + 0.025*x**2 - 5.323*x**3) if x > 0 else (1.362 - 0.14*x)

	lrhk = np.mean(system.extras.rhk)
	y = 5. + lrhk
	## Noyes (1984), Eq 3
	log_P = (0.324 - 0.4*y - 0.283*y**2 - 1.325*y**3) + log_tau

	msg = yellow('RESULT: ') + 'from Noyes (1984), Prot = %f d' % (10**log_P)
	clogger.info(msg)


	# ## Mamajek & Hillenbrand (2008), Eq 3
	# log_tau = -38.053 - 17.912*lrhk - 1.6675*lrhk**2
	# tau = 10**log_tau

	# ## Mamajek & Hillenbrand (2008), Eq 5
	# P = (0.808 - 2.966*(lrhk+4.52)) * tau

	# msg = yellow('RESULT: ') + 'from Mamajek & Hillenbrand (2008), Prot = %f d' % (P)
	# clogger.info(msg)



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


def do_clean(system):
	""" Run the CLEAN algorithm """

	if not periodogram_CLEAN_available:
		msg = red('ERROR: ') + 'This extension is not available. Something went wrong on install...'
		clogger.fatal(msg)
		return

	time, rv, err = system.time, system.vrad, system.error
	plow = 1.0

	msg = blue('INFO: ') + 'Running CLEAN...'
	clogger.info(msg)

	# run CLEAN
	df, Nf, niter, D_r, D_i, W_r, W_i, C_r, C_i = clean(time, rv, err, plow)

	# get the amplitudes of the complex CLEANed spectrum
	C = C_r * C_r + C_i * C_i

	msg = yellow('    : ') + 'number of iterations: %d\n' % niter
	msg += yellow('    : ') + 'CLEAN spectrum resolution: %8.3e [day^-1]\n' % df
	msg += yellow('    : ') + 'CLEAN spectrum Pmin: %5.3f [days]\n' % plow
	msg += yellow('    : ') + 'Normalizing CLEAN spectrum by %8.3e\n' % (max(C),)
	clogger.info(msg)

	# Define the frequency ranges for Dnu/Cnu and Wnu
	nu1 = np.arange(-Nf*df, Nf*df+df, df)
	# nu2 = np.arange(-Nf*df*2, 2*Nf*df+df, df)

	if len(nu1) == len(C) + 1:
		nu1 = nu1[:-1]

	# (arbitrary) renormalization
	C = C / max(C)

	plt.figure()
	ax1 = plt.subplot(111)
	ax1.semilogx((1./nu1), C,'k-')
	ax1.set_xlabel('Period [d]')
	ax1.set_ylabel('Arbitrary power')
	plt.show()



def do_create_planets(s):
	msg = blue('INFO: ') + 'Starting the planet creator'
	clogger.info(msg)

	if s != '':  # parse options from the command line
		import re
		s = s.lower()
		try:
			# match periods
			regex1 = re.compile("p\((.*?)\)")
			p = regex1.findall(s)[0]
			periods = [float(i) for i in p.split(',')]
			# match eccentricities
			regex2 = re.compile("e\((.*?)\)")
			e = regex2.findall(s)[0]
			eccentricities = [float(i) for i in e.split(',')]
			# match semi-amplitudes
			regex3 = re.compile("k\((.*?)\)")
			k = regex3.findall(s)[0]
			semi_amplitudes = [float(i) for i in k.split(',')]
		except IndexError:
			msg = red('ERROR: ') + "Couldn't parse arguments"
			clogger.fatal(msg)
			return

		if not len(periods) == len(eccentricities) == len(semi_amplitudes):
			msg = red('ERROR: ') + 'Non-matching number of options'
			clogger.fatal(msg)
			return

		nplanets = len(periods)

		try:
			regex = re.compile("n\((.*?)\)")
			nobs = int(regex.findall(s)[0])
		except IndexError:
			nobs = np.random.randint(30, 220)

		try:
			regex = re.compile("file\((.*?)\)")
			save_filename = regex.findall(s)[0]
		except IndexError:
			save_filename = None

		type_noise = 'white'  # for now...
		filename = None

	else:

		try:
			nplanets = int(raw_input('How many planets? '))
		except ValueError:
			msg = red('ERROR: ') + "I don't know how many planets to make\n"
			clogger.fatal(msg)
			msg = blue('INFO: ') + 'Finished planet creator\n'
			clogger.info(msg)
			return

		if ask_yes_no('Specify the periods (y/N)? ', default=False):
			periods = []
			for i in range(nplanets):
				p = float(raw_input('\tperiod %d:  ' % (i+1)))
				periods.append(p)
		else:
			periods = None

		if ask_yes_no('Specify the eccentricities (y/N)? ', default=False):
			eccentricities = []
			for i in range(nplanets):
				e = float(raw_input('\teccentricity %d:  ' % (i+1)))
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
			nobs = np.random.randint(30, 220) if obs in ('r', '') else int(obs)

		noise = 'wrong_option'
		while noise not in ('no', 'w', '', 'r', 'R', 'wr', 'rw'):
			noise = raw_input('What type of noise? (no / W / r / wr) ')
			if noise in ('w', ''):
				type_noise = 'white'
			elif noise in ('r', 'R'):
				type_noise = 'correlated'
			elif noise in ('wr', 'rw'):
				type_noise = 'white + correlated'
			elif noise in ('no'):
				type_noise = 'no'
			else:
				print "I don't understand that. Try 'no', w', 'r' or 'wr'."

	print
	msg = blue('INFO: ') + 'Generating %d planet(s)\n' % (nplanets)
	if filename is not None:
		msg += blue('    : ') + '-> sampling from %s\n' % (filename)
	else:
		msg += blue('    : ') + '-> randomly spaced sampling\n'
	msg += blue('    : ') + '-> %d observations\n' % (nobs)
	msg += blue('    : ') + '-> %s noise\n' % (type_noise)
	clogger.info(msg)

	# get rv curve with random parameters
	def gen_rv(nplanets, nobs, sampling_file, periods, eccentricities, type_noise, temp):

		if sampling_file is not None:  # the sampling is set in the file
			times_sampled_from_file = np.loadtxt(sampling_file, usecols=(0,), skiprows=2)
			times = times_sampled_from_file
			times_full = np.linspace(min(times_sampled_from_file), max(times_sampled_from_file), 1000)
		else:  # the sampling is random
			times_full = np.linspace(2449460, 2452860, 1000)
			# sample (N points) randomly from the data to produce unevenly sampled time series
			rand_args = [i for i in sorted(random.sample(xrange(len(times_full)), nobs))]
			times = times_full[rand_args]

		vel_full = zeros_like(times_full)  # rv at oversampled times
		vel_full_each = zeros_like(times_full)  # rv at oversampled times for each planet
		vel_all = np.zeros((len(times), nplanets))  # rv at sample times for each planet
		vel_each = zeros_like(times)  # rv at sample times for each planet (temporary)
		vel_total = zeros_like(times)  # rv at sample times

		output = '# planet %d with P=%f, K=%f, e=%f, w=%f, t0=%f\n'
		for planet in range(nplanets):
			if periods is None: 
				P = np.random.rand()*998. + 2.  # random period, U(2, 998)
			else: 
				P = periods[planet]  # user-provided period

			if eccentricities is None:
				ecc = np.random.beta(0.867, 3.03)  # from Kipping 2013
			else:
				ecc = eccentricities[planet]  # user-provided eccentricity

			K = np.random.rand()*10. + 1.

			omega = np.random.rand()*2.*pi 

			if type_noise is 'white':
				noise = np.random.randn(len(vel_each))  # sigma??????
			else:
				noise = np.ones_like(vel_each)

			# log
			temp.write(output % (planet+1, P, K, ecc, omega, 2452000))

			# get RVs for this planet
			get_rvn(times, P, K, ecc, omega, 2452000, 0., vel_each)
			if type_noise is 'white':
				vel_each = vel_each + noise  # add noise

			# get RVs for this planet at oversampled times
			get_rvn(times_full, P, K, ecc, omega, 2452000, 0., vel_full_each)

			# store
			vel_all[:,planet] = vel_each
			vel_total += vel_each
			vel_full += vel_full_each

		return times_full, vel_full, times, vel_total, vel_all, noise

	# wrapper function that takes care of adding noise and logging to file
	def generate_planets(nplanets, nobs, sampling_file, saving_file, type_noise, periods, eccentricities):
		if saving_file is None:
			# create temporary file to store simulation info
			tf = tempfile.NamedTemporaryFile(dir='./')
		else:
			tf = open(saving_file, 'w')

		tf.write('## Planets created on ')
		tf.write(time.strftime("%d/%m/%Y - %H:%M:%S\n#\n"))

		times_full, vel_full, times, vel, vel_all, noise = \
			gen_rv(nplanets, nobs, sampling_file, periods, eccentricities, type_noise, tf)

		tf.write('#\n## %d observations, %s noise\n#\n' % (nobs, type_noise))
		tf.write('# bjd \t vrad(km/s) \t svrad(km/s)\n')

		for t, v, e in zip(times, vel, noise):
			tf.write('%f\t%f\t%f\n' % (t, v*1e-3, abs(e*1e-3)))

		try:
			tf.delete = False
		except AttributeError:
			pass

		print 'Output to %s' % tf.name
		tf.close()

		return times_full, vel_full, times, vel, vel_all

	save_filename = None
	r = generate_planets(nplanets, nobs, filename, save_filename, type_noise, periods, eccentricities)
	times_full, vel_full, times, vel1, vel2 = r


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
