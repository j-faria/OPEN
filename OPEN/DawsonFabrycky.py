# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

import sys,os,string,dircache
from pylab import *
from math import pi
from math import exp
from math import sqrt
from math import atan2
from math import cos
from math import sin
import cmath
import random
from ext import periodogram_DF


ofac = 2.0
plow = 0.9

def do_Dawson_Fabrycky(system):
	time, rv, err = system.time, system.vrad, system.error
	# Make permutations for indices of length n

	def permutate(n):

		A=[]
		for i in range(n):
			a = int(n*random.random())
			if a==n: a=n-1
			A.append(a)

		return A


	# Calculate the spectral window function and its phase angles

	def specwindow(freq,time):

		n = len(time)
		W = [sum([cmath.exp(-2.j*pi*f*t) for t in time])/float(n) for f in freq]

		amp = [sqrt(t.real*t.real + t.imag*t.imag) for t in W]
		phase = [atan2(t.imag,t.real) for t in W]

		return amp,phase


	def fitcirc(per,K,Tc,gam,t):

		phase = 2*pi*Tc/per

		omega = 2. * pi / per

		rv = K * cos(omega * t - phase) + gam

		return rv



	def fitkep(per,ecc,pomega,K,T0,gam,t):


		afit = K*cos(pomega/180*pi)
		bfit = -K*sin(pomega/180*pi)
		cfit = K*ecc*cos(pomega/180*pi) + gam

		M=(t-T0)/per*2.*pi

		sE = 0.
		sEtemp = 2.
		while abs(sE - sEtemp) > 0.05:

			sEtemp = sE
			EA = ecc * sE + M
			sE = sin(EA)
			
		cE = cos(EA)

		cosnu = (cE - ecc) / (1. - ecc * cE)
		sinnu = sqrt(1. - ecc * ecc) * sE / (1. - ecc * cE)

		rv = afit * cosnu + bfit * sinnu + cfit

		return rv

	plow = 0.5

	n = len(err)

	### GET THE REAL PERIODOGRAM
	freq,power,a_cos,b_sin,c_cte,phi,fNy,xdif = periodogram_DF.periodogram_DF(array(time),array(rv),array(err),ofac,plow)

	### GET THE WINDOW FUNCTION AT THOSE FREQUENCIES + plot dials
	amp,phase = specwindow(freq,time)

	figure(num=1,figsize=(11.69,8.27))


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

	rv_fake1 = [a1*cos(fmax1*2.*pi*i) + b1*sin(fmax1*2.*pi*i) + c1 for i in timefake]
	rv_fake2 = [a2*cos(fmax2*2.*pi*i) + b2*sin(fmax2*2.*pi*i) + c2 for i in timefake]
	rv_fake3 = [a3*cos(fmax3*2.*pi*i) + b3*sin(fmax3*2.*pi*i) + c3 for i in timefake]
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

	freq,power,a_cos,b_sin,c_cte,phi,fNy,xdif = periodogram_DF.periodogram_DF(array(timefake),array(rv_fake1),array(errfake),ofac,plow)

	#nf = len(freq)
	#temp = sorted(zip(power,freq,a_cos,b_sin,c_cte,phi))
	#fmax1 = (temp[nf-1])[1]

	#fmax2 = fmax1
	#i = 1
	#while abs(fmax2 - fmax1) < 0.01:
		
	#	i += 1
	#	fmax2 = (temp[nf-i])[1]

	#fmax3 = fmax2

	#j = i
	#while abs(fmax3 - fmax2) < 0.01 or abs(fmax3 - fmax1) < 0.01:
		
	#	j += 1
	#	fmax3 = (temp[nf-j])[1]


	#print 'Maximum frequencies at: ',fmax1,fmax2,fmax3
	#print 'That is periods at: ',round(1./fmax1,2),round(1./fmax2,2),round(1./fmax3,2)

	#ph1 = (temp[nf-1])[5]
	#ph2 = (temp[nf-i])[5]
	#ph3 = (temp[nf-j])[5]

	ind1 = list(freq).index(fmax1)
	ind2 = list(freq).index(fmax2)
	ind3 = list(freq).index(fmax3)

	subplot(4,1,2)

	semilogx([1./i for i in freq],power,'k-')

	semilogx(1./fmax1,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax2,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
	semilogx(1./fmax3,max(power)+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)

	semilogx([1./fmax1,1./fmax1+0.045*cos(phi[ind1])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind1])],'k-',lw=1)
	semilogx([1./fmax2,1./fmax2+0.045*cos(phi[ind2])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind2])],'k-',lw=1)
	semilogx([1./fmax3,1./fmax3+0.045*cos(phi[ind3])],[max(power)+0.1,max(power)+0.1+0.045*sin(phi[ind3])],'k-',lw=1)


	xlim(plow,xdiff*ofac)
	ylim(0.0,max(power)+0.2)


	#### 2nd FAKE

	freq,power,a_cos,b_sin,c_cte,phi,fNy,xdif = periodogram_DF.periodogram_DF(array(timefake),array(rv_fake2),array(errfake),ofac,plow)


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

	freq,power,a_cos,b_sin,c_cte,phi,fNy,xdif = periodogram_DF.periodogram_DF(array(timefake),array(rv_fake3),array(errfake),ofac,plow)


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
	# close()





