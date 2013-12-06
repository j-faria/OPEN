# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
from numpy import polyfit, RankWarning, append
# see http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
from numpy import zeros_like, savetxt
from matplotlib.pylab import *
import warnings
from datetime import datetime
import subprocess
from logger import clogger, logging
import sys, os
from scipy.optimize import leastsq
from ext.get_rv import get_rv
from galileo import *
from shell_colors import yellow, red, blue

def do_fit(system, verbose):
	try:
		degree = system.model['d']
	except TypeError:
		msg = red('Error: ') + 'Need to run mod before fit. '
		clogger.error(msg)
		return
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


def do_genetic(system):
	vel = zeros_like(system.time)

	def chi2_1(params):
		P = params[0]
		K = params[1]
		ecc = params[2]
		omega = params[3]
		t0 = params[4]
		get_rv(system.time, P, K, ecc, omega, t0, vel)
		return 1./sum(((system.vrad - vel)/system.error)**2)

    #create an initial population of 10 chromosomes
	p = Population(10)
	#use fitness (above) as our evaluation function
	p.evalFunc = chi2_1
	#minimum values the genes can take
	p.chromoMinValues = [5, 0, 0, 0, 2451000]    
	#maximum values the genes can take
	p.chromoMaxValues = [10000, 100, 0.9, 6.28, 2451500]
	#use integers instead of floats
	p.useInteger = 0
	#always crossover
	p.crossoverRate = 1.0
	#mutate, but not very often
	p.mutationRate = 0.05
	#use roulette (monte carlo) selection
	p.selectFunc = p.select_Roulette
	#use a full replacement size
	p.replacementSize = p.numChromosomes
	#use one point crossover
	p.crossoverFunc = p.crossover_OnePoint
	#use the default mutation routines
	p.mutateFunc = p.mutate_Default
	#use steady-state replacement with no duplication
	p.replaceFunc = p.replace_SteadyStateNoDuplicates
	#p.replaceFunc = p.replace_Generational

	#finish initializing the population. THIS MUST BE CALLED after settings the
	#variables above, but before actually running the GA!
	p.prepPopulation()

	for i in range(500):
	  #evaluate each chromosomes
	  p.evaluate()
	  #apply selection
	  p.select()
	  #apply crossover
	  p.crossover()
	  #apply mutation
	  p.mutate()
	  #apply replacement
	  p.replace()
	  #print info
	#  print p.minFitness, p.maxFitness, p.avgFitness, p.sumFitness
	print 'Genetic:', p.bestFitIndividual, p.bestFitIndividual.fitness
	lm = do_lm(system, p.bestFitIndividual.genes)
	lm_par = lm[0]
	print 'LM:', lm_par

	# get best solution curve
	new_time = system.get_time_to_plot()
	vel = zeros_like(new_time)

	P, K, ecc, omega, t0 = p.bestFitIndividual.genes
	get_rv(new_time, P, K, ecc, omega, t0, vel)
	# plot RV with time
	plot(system.time, system.vrad, 'o')
	plot(new_time, vel, '-')

	P, K, ecc, omega, t0 = lm_par
	get_rv(new_time, P, K, ecc, omega, t0, vel)
	plot(new_time, vel, 'r-')
	show()

	return

def do_lm(system, x0):

	vel = zeros_like(system.time)

	def chi2_2(params):
		P = params[0]
		K = params[1]
		ecc = params[2]
		omega = params[3]
		t0 = params[4]
		get_rv(system.time, P, K, ecc, omega, t0, vel)
		return system.vrad - vel

	return leastsq(chi2_2, x0, full_output=0)


def do_multinest(system):
	msg = blue('INFO: ') + 'Transfering data to MultiNest...'
	clogger.info(msg)

	# write data to file to be read by MultiNest
	nest_filename = 'input.rv'
	nest_header = 'file automatically generated for MultiNest analysis, ' + str(datetime.now())
	savetxt(nest_filename, zip(system.time, system.vrad, system.error),
		    #header=nest_header,
		    fmt=['%12.4f', '%6.1f', '%6.2f'])


	msg = blue('INFO: ') + 'Starting MultiNest...'
	clogger.info(msg)

	cmd = 'mpirun -np 2 ./OPEN/multinest/gaussian'
	subprocess.call(cmd, shell=True)
	# os.system(cmd)


	return