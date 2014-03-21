# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

# standard library imports
import os
import glob
from fileinput import FileInput
from itertools import islice, chain
from collections import namedtuple # this requires Python >= 2.6
# other imports
import numpy as np
import matplotlib.pylab as plt
import scipy.stats
from scipy.signal import argrelmax
# from matplotlib import rc, rc_params_from_file
import matplotlib.gridspec as gridspec

# intra-package imports
try:
    import pyqtgraph as pg
except Exception:
    print 'pyqtgraph will not be installed'

import rvIO
from .utils import unique
from .logger import clogger, logging
from shell_colors import yellow, blue
from .utils import day2year, rms, ask_yes_no, triangle_plot
from ext.get_rvN import get_rvn


# if the interpreter calls open.py from another directory, matplotlib will not
# access the matplotlibrc file and will revert to its defaults. This makes sure
# we use it!
# p = os.path.dirname(os.path.realpath(__file__)) # path to this file (classes.py)
# rc( rc_params_from_file( p+'/../matplotlibrc' ) )
plt.ion()

class rvSeries:
    """
    A container class that holds the observed radial velocity data. 
    
    It includes a *plot* method that takes care of points from one or more 
    files.
            
    Parameters
    ----------
    filenames: string
        One or more files to read radial velocities from.
    """

    def __repr__(self):
        return "System with RV info from %s" % self.provenance
    # def __str__(self):
    #     return "member of Test"

    def __init__(self, *filenames, **kwargs):
        from shell_colors import red
        try:
            assert len(filenames)>=1
        except AssertionError:
          msg = red('ERROR: ') + 'Need at least one file to read! (aka, does that file exist?)'
          clogger.fatal(msg)
          raise IOError

        # don't repeat files
        filenames = unique(filenames)
          
        # skip header lines? (no by default)
        try:
          skip = kwargs.pop('skip')
        except KeyError:
          skip = 0


        # read data
        try:
          data, self.provenance = rvIO.read_rv(*filenames, skip=skip)
        except ValueError as e:
          print e
          msg = red('ERROR: ') + 'If your files have header lines set --skip option.\n'
          clogger.fatal(msg)
          return

        
        self.time, self.vrad, self.error = data.pop('jdb'), data.pop('vrad'), data.pop('svrad')

        # save the extra quantities as a namedtuple if we read them
        extra = namedtuple('Extra', data.keys(), verbose=False)
        self.extras = extra(**data)
        self.extras_names = data.keys()

        # time, vrad and error can change, 
        # the *_full ones correspond always to the full set
        self.time_full = self.time
        self.vrad_full = self.vrad
        self.error_full = self.error
        # same for extras
        self.extras_full = self.extras

    # associated model to be adjusted to the data
    # this will be a dictionary with the following key:value pairs:
    #   k : number of keplerians in the model
    #   d : degree of drift in the model
    #   drift : array with parameters from drift fit (the polynomial coefficients)
    model = None

    # associated fit to be carried out with the genetic / LM algorithms
    # this will be a dictionary with the following key:value pairs:
    #   params : final parameters (5*model['k']+1)
    #   residuals : residuals of the fit
    #   chi2 : reduced(!) chi square value of the fit
    fit = None

    def save(self, filename):
        rvIO.write_rv(self, filename)

    def stats(self):
        # verbose stats about data
        t, rv = self.time, self.vrad
        info = blue('INFO: ') 
        sinfo = blue('    : ') 
        stats = None
        tspan = max(t) - min(t)
        rvspan = max(rv) - min(rv)
        stats = '\n'
        stats += info + "Timespan : %f days = %f years   ---   %f JD, %f JD\n" % (tspan, day2year(tspan), min(t), max(t))
        stats += sinfo + "RV span  : %f km/s = %f m/s\n" % (rvspan, rvspan*1e3)
        stats += sinfo + "RV rms [m/s] : %f\n\n" % rms(rv*1e3)
        stats += sinfo + "{:14s} : {:10.3f}\n".format('<RV> [m/s]', np.mean(rv*1e3))
        if len(self.extras) == 0: 
            clogger.info(stats)
            return

        stats += sinfo + "{:14s} : {:10.3f}\n".format('<fwhm> [km/s]', np.mean(self.extras.fwhm))
        stats += sinfo + "{:14s} : {:10.3f}\n".format('<contrast>', np.mean(self.extras.contrast) )
        stats += sinfo + "{:14s} : {:10.3f}\n".format('<BIS> [km/s]', np.mean(self.extras.bis_span) )
        stats += sinfo + "{:14s} : {:10.3f}\n".format('<S_index> [MW]', np.mean(self.extras.s_mw) )
        stats += sinfo + "{:14s} : {:10.3f}\n".format('<log(rhk)>', np.mean(self.extras.rhk))

        clogger.info(stats)

    def do_plot_obs(self, newFig=True, leg=True):
        """ Plot the observed radial velocities as a function of time.
        Data from each file are color coded and labeled.
        """

        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = self.time, self.vrad, self.error # temporaries
        
        if newFig: 
            plt.figure()
        # plot each files' values
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
            #                     height=err[:m], beam=0.5,\
            #                     pen=pg.mkPen(None))
                                # pen={'color': 0.8, 'width': 2})
            # p.addItem(e)
            # p.plot(t[:m], rv[:m], symbol='o')
            plt.errorbar(t[:m], rv[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err = t[m:], rv[m:], err[m:]
        
        plt.xlabel('Time [days]')
        plt.ylabel('RV [km/s]')
        if leg: plt.legend()
        plt.tight_layout()
        plt.ticklabel_format(useOffset=False)
        plt.show()
        # pg.QtGui.QApplication.exec_()

    def do_plot_drift(self):
        """ Plot the observed radial velocities as a function of time, plus an
        extra drift of specified degree (see *mod*). Lower panel presents RV 
        minus drift. Data from each file are color coded and labeled.
        """
        try:
            drift = self.model['drift']
        except KeyError: # there is no drift set
            return 

        colors = 'rgbmk' # possible colors
        t, rv, err = self.time, self.vrad_full, self.error # temporaries
        
        plt.figure()

        ax1 = plt.subplot(2,1,1)
        # plot each file's values
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            ax1.errorbar(t[:m], rv[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err = t[m:], rv[m:], err[m:]
        
        p = numpy.polyval(drift, self.time) # normal polynomial, for 2nd plot
        st = np.sort(self.time) # need this otherwise the plot gets messed up
        sp = np.polyval(drift, st) # "sorted" polynomial
        ax1.plot(st, sp, 'y-')

        t, rv, err = self.time, self.vrad, self.error # temporaries
        ax2 = plt.subplot(2,1,2, sharex=ax1)
        # plot each file's values minus the drift
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            ax2.errorbar(t[:m], rv[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err, p = t[m:], rv[m:], err[m:], p[m:]

        plt.xlabel('Time [days]')
        plt.ylabel('RV [km/s]')
        # plt.legend()
        plt.tight_layout()
        plt.ticklabel_format(useOffset=False)
        plt.show()

    def do_plot_extras(self, extra):
        """ Plot other observed quantities as a function of time.

        Parameterslistcommands
        ----------
        extra: string
          One of the quantities available in system.extras
        """
        # import pyqtgraph as pg

        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t = self.time

        # handle inexistent field
        if extra not in self.extras._fields:
          from shell_colors import red
          msg = red('ERROR: ') + 'The name "%s" is not available in extras.\n' % extra
          clogger.fatal(msg)
          return

        i = self.extras._fields.index(extra) # index corresponding to this quantity

        plt.figure()
        # p = pg.plot()
        if extra == 'rhk':
            plt.errorbar(t, self.extras.rhk, yerr=self.extras.sig_rhk, fmt='o', label=extra)
        else:
            plt.plot(t, self.extras[i], 'o', label=extra)
        
        plt.xlabel('Time [days]')
        plt.ylabel(extra + ' []')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05))
        plt.minorticks_on()
        plt.tight_layout()
        plt.show()
        # pg.QtGui.QApplication.exec_()

    def do_plot_fit(self):
        """ Plot the observed radial velocities together with the 
        current best fit curve as well as phased RV curves for each 
        planet considered in the fit. Data from each file are color 
        coded and labeled. 
        """
        if self.fit is None: return # there is no fit yet

        colors = 'rgbmk' # possible colors
        t, rv, err = self.time, self.vrad, self.error # temporaries

        tt = np.linspace(self.time.min()-200, self.time.max()+300, 400)
        final = np.zeros_like(tt)
        
        par = self.fit['params']
        keplerians = int(len(par)/6)
        P, K, ecc, omega, T0, gam = [par[j::6] for j in range(6)]
        gam = gam[0]

        get_rvn(tt, P, K, ecc, omega, T0, gam, final)

        ### observations + fitted curve + residuals
        plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
        ax1 = plt.subplot(gs[0])
        # plot each file's values
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            ax1.errorbar(t[:m], rv[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err = t[m:], rv[m:], err[m:]
        ax1.plot(tt, final, 'k-')
        ax1.ticklabel_format(useOffset=False)

        # redo this...
        t, rv, err = self.time, self.vrad, self.error # temporaries
        final = np.zeros_like(t)
        get_rvn(t, P, K, ecc, omega, T0, gam, final)

        ax2 = plt.subplot(gs[1], sharex=ax1)
        # plot residuals
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            ax2.errorbar(t[:m], rv[:m]-final[:m], yerr=err[:m], \
                       fmt='o'+colors[i], label=fname)
            t, rv, err, final = t[m:], rv[m:], err[m:], final[m:]
        ax2.axhline(y=0, xmin=0, xmax=1, ls='--', color='k')

        # redo this...
        t, rv, err = self.time-self.time.min(), self.vrad, self.error # temporaries
        tt = np.linspace(self.time.min(), self.time.max(), 400)
        final = np.zeros_like(tt)

        ### phased RV curves
        plt.figure()
        gs = gridspec.GridSpec(keplerians, 1)
        for i in range(keplerians):
          get_rvn(tt, P[i], K[i], ecc[i], omega[i], T0[i], gam, final)
          ax = plt.subplot(gs[i])
          # plot each file's values
          for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
              m = n-nout # how many values are there after restriction
              ax.errorbar(np.modf(abs(t[:m]-T0[i])/P[i])[0], rv[:m], yerr=err[:m], \
                           fmt='o'+colors[i], label=fname)
              t, rv, err = t[m:], rv[m:], err[m:]
          ax.plot( np.modf(abs(tt-T0[i])/P[i])[0], final, 'k-')

        plt.show()

    def save_fit(self, params, chi2):
        """ Helper function to save results from a fit to the system """
        if self.fit is not None: # there is already a fit
          clogger.warning(yellow('Warning: ')+'Replacing older fit')

        self.fit = {}
        self.fit['params'] = params
        self.fit['chi2'] = chi2

        P, K, ecc, omega, T0, gam = [params[i::6] for i in range(6)]
        gam = gam[0]

        final = np.zeros_like(self.time)
        get_rvn(self.time, P, K, ecc, omega, T0, gam, final)
        self.fit['residuals'] = self.vrad - final

    def get_nyquist(self, smallest=False):
        """
        Calculate the average Nyquist frequency. 
        If smallest is True, calculates highest Nyquist frequency 
        (with smallest time interval)
          
        Returns
        -------
        Nyquist frequency : float
            Half the sampling frequency of the time series.
        """
        # this way is about 2 times faster than np.ediff1d
        if smallest:
            return 0.5*1./np.min(self.time[1::]-self.time[0:-1])
        else:
            return 0.5*1./np.mean(self.time[1::]-self.time[0:-1])

    def get_time_to_plot(self):
        """
        Returns sufficiently resolved time vector for plots
        """
        std = self.time.std()
        N = len(self.time)
        minim = self.time.min() - 1.*std
        maxim = self.time.max() + 1.*std
        return np.linspace(minim, maxim, 10*N)
        
        
class PeriodogramBase:
    """
    Base class for all periodograms.
    
    It contains a standard *plot* method to visualize periodograms and methods 
    to calculate significance, FAPs and other periodogram-related statistics.

    Attributes
    ----------
    power : array
        The periodogram power.
    freq : array
        The frequencies at which the power are evaluated.
    """
    
    name = None # which periodogram


    def get_peaks(self, n=1, output_period=False):
        """ Get the frequencies of the largest n peaks in the periodogram. 
        If output_period is True, also return the periods. 
        """

        nf = len(self.freq)
        if self.name is 'SpectralWindow':
            self.power = self.amp

        ind = argrelmax(self.power, order=20)[0]
        indsort = np.argsort(self.power[ind])

        peaks = np.array([self.freq[ind[indsort[-i]]] for i in range(1,n+1)])
        if output_period:
            return [(1/peak, peak) for peak in peaks]
        else:
            return peaks

        # temp = sorted(zip(self.power, self.freq))

        # # largest peak
        # fmax1 = (temp[nf-1])[1]
        # if n==1:
        #     if output_period: 
        #         return fmax1, 1./fmax1
        #     else:
        #         return fmax1

        # # second largest
        # fmax2 = fmax1
        # i = 1
        # while abs(fmax2 - fmax1) < 0.01:
        #     i += 1
        #     fmax2 = (temp[nf-i])[1]

        # if n==2:
        #     if output_period: 
        #         return [(fmax1, 1./fmax1), 
        #                 (fmax2, 1./fmax2)]
        #     else:
        #         return fmax1, fmax2

        # # third largest
        # fmax3 = fmax2
        # j = i
        # while abs(fmax3 - fmax2) < 0.01 or abs(fmax3 - fmax1) < 0.01:
        #     j += 1
        #     fmax3 = (temp[nf-j])[1]

        # if n==3:
        #     if output_period: 
        #         return [(fmax1, 1./fmax1), 
        #                 (fmax2, 1./fmax2), 
        #                 (fmax3, 1./fmax3)]
        #     else:
        #         return fmax1, fmax2, fmax3

    def FAP(self, Pn):
      """
      Calculate the false-alarm probability (FAP).

      The FAP denotes the probability that at least one out of M
      independent power values in a prescribed search band of a
      power spectrum computed from a white-noise time series is
      as large as or larger than the threshold Pn.
      It is assessed through
      
      .. math:: FAP(Pn) = 1 - (1-Prob(P>Pn))^M \\; ,
      
      where "Prob(P>Pn)" depends on the type of periodogram and 
      normalization and is calculated by using the *prob* method;
      M is the number of independent power values and (an estimate)
      is computed internally (see self.M).

      Parameters
      ----------
        Pn : float
            Power threshold.
        
      Returns
      -------
        FAP : float
            False alarm probability.
      """
      prob = self.M * self.prob(Pn)
      if prob > 0.01:  return 1.-(1.-self.prob(Pn))**self.M
      return prob

    def powerLevel(self, FAPlevel):
      """
        Power threshold for FAP level.
      
        Parameters
        ----------
        FAPlevel : float or array
              "False Alarm Probability" threshold
      
        Returns
        -------
        Threshold : float or array
            The power threshold pertaining to a specified
            false-alarm probability (FAP). Powers exceeding this
            threshold have FAPs smaller than FAPlevel.
      """
      Prob = 1.-(1.-FAPlevel)**(1./self.M)
      return self.probInv(Prob)  

    def FAP_by_bootstrap(self):

        name = '_' + self.__class__.__name__
        # access the instance's __calcPeriodogramFast function
        exec 'calc = self.' + name + '__calcPeriodogramFast'

        # temporaries
        f = self.freq
        p = self.power

        perc01 = 0.001 # 0.1% FAP
        perc1 = 1.  # 1% FAP
        perm = int(1000/perc1) # (use 1000 for 1% fap or 10000 for 0.1% fap)

        maxPowers = []
        for k in tqdm(range(perm)):
            permutted = np.random.permutation(zip(self.y, self.error))
            self.y = permutted[:,0]
            self.error = permutted[:,1]

            calc()
            # self._plot_pg()
            powermaxP = self.power.max()
            maxPowers.append(powermaxP)
        #print maxPowers

        peaks = np.sort(maxPowers)
        # index01 = int( ((1-perc01/100.0) * len(peaks)) )
        index1 = int( ((1-perc1/100.0) * len(peaks)) )
        # powerFAP_01 = peaks[index01]
        powerFAP_1 = peaks[index1]


        plt.semilogx(1./f, p, 'k-')
        # plt.axhline(powerFAP_01,c='r',ls=':')
        plt.axhline(powerFAP_1,c='r',ls='--')
        plt.show()
        #         if orbit == 'circ':
        #             powermaxP = (periodogram.periodogram(bjd,data_perm,sigma_perm,ofac,plow))[3]
        #         if orbit == 'kep':
        #             powermaxP = (periodogram_kep.periodogram_kep(bjd,data_perm,sigma_perm,ofac,plow))[3]
        # #       print k
        #         maxPowers.append(powermaxP)

    def _plot(self, doFAP=False, faps=None, verts=None, newFig=True, axes=None, **kwargs):
      """
        Create a plot.
      """
      xlabel = 'Period [d]'
      ylabel = 'Power'

      if newFig: 
        self.fig = plt.figure()

      if axes is None:
        self.ax = self.fig.add_subplot(1,1,1)
      else:
        self.ax = axes

      self.ax.set_title("Normalized periodogram")
      self.ax.set_xlabel(xlabel)
      self.ax.set_ylabel(ylabel)
      if self.power.max() < 1e-300:  # apparently, Metplotlib can't plot these small values
        clogger.warning(yellow('Warning: ')+'Max value < 1e-300, plotting normalized periodogram')
        self.ax.semilogx(1./self.freq, self.power/self.power.max(), 'b-', **kwargs)
      else:
        self.ax.semilogx(1./self.freq, self.power, 'b-', **kwargs)
      # plot FAPs
      if doFAP and (faps is None):
        clogger.warning(yellow('Warning: ')+'Plotting default FAPs')
      if doFAP: # do default FAPs of 10%, 1% and 0.1%
        pmin = 1./self.freq.min()
        pmax = 1./self.freq.max()
        plvl1 = self.powerLevel(0.1) # 10% FAP
        plvl2 = self.powerLevel(0.01) # 1% FAP
        plvl3 = self.powerLevel(0.001) # 0.1% FAP
        self.ax.axhline(y=plvl1, color='k', ls='-', label='10%')
        self.ax.axhline(y=plvl2, color='k', ls='--', label='1%')
        self.ax.axhline(y=plvl3, color='k', ls=':', label='0.1%')
        self.ax.legend(frameon=True)

      # plot vertical lines
      if verts is not None:
        for v in verts:
          self.ax.axvline(x=v, color='k', ls='--', lw=2, alpha=0.5) 
          # if v==18:
          #   self.ax.axvline(x=v, color='r', ls='--', lw=2) 


      plt.tight_layout()

    def _plot_pg(self, doFAP=False, verts=None, newFig=True, axes=None):

      if self.name == 'Bayesian Lomb-Scargle':
        p = pg.plot(1./self.freq, self.power/self.power.max(), title="Periodogram")
      else:
        p = pg.plot(1./self.freq, self.power, title="Periodogram")
      p.plotItem.setLogMode(True, False)
      pg.QtGui.QApplication.exec_()

    def _output(self, verbose=False):
      """
      Some statistical output.
      """
      from shell_colors import blue
      # Index with maximum power
      bbin = argmax(self.power)
      # Maximum power
      pmax = self._upow[bbin]

      rms = sqrt(self._YY * (1.-pmax))

      # Get the curvature in the power peak by fitting a parabola y=aa*x^2
      if (bbin > 1) and (bbin < len(self.freq)-2):
        # Shift the parabola origin to power peak
        xh = (self.freq[bbin-1:bbin+2] - self.freq[bbin])**2
        yh = self._upow[bbin-1:bbin+2] - self._upow[bbin]
        # Calculate the curvature (final equation from least square)
        aa = sum(yh*xh)/sum(xh*xh)
        nt = float(self.N)
        f_err = sqrt(-2./nt * pmax/aa*(1.-pmax)/pmax)
        Psin_err = sqrt(-2./nt* pmax/aa*(1.-pmax)/pmax) / self.freq[bbin]**2
      else:
        f_err = None
        Psin_err = None

      fbest = self.freq[bbin]
      amp = sqrt(self._a[bbin]**2 + self._b[bbin]**2)
      ph  = arctan2(self._a[bbin], self._b[bbin]) / (2.*pi)
      T0  = min(self.th) - ph/fbest
      # Re-add the mean
      offset = self._off[bbin] + self._Y

      # Statistics
      print "Generalized LS - statistical output"
      print 33*"-"
      if verbose:
        print "Number of input points:     %6d" % (nt)
        print "Weighted mean of dataset:   % e" % (self._Y)
        print "Weighted rms of dataset:    % e" % (sqrt(self._YY))
        print "Time base:                  % e" % (max(self.th) - min(self.th))
        print "Number of frequency points: %6d" % (len(self.freq))
        print
        print "Maximum power, p :    % e " % (self.power[bbin])
        print "Maximum power (without normalization):   %e" % (pmax)
        print "Normalization    : ", self.norm
        print "RMS of residuals :    % e " % (rms)
        if self.error is not None:
          print "  Mean weighted internal error:  % e" %(sqrt(nt/sum(1./self.error**2)))
        print "Best sine frequency : % e +/- % e" % (fbest, f_err)
        print "Best sine period    : % e +/- % e" % (1./fbest, Psin_err)
        print "Amplitude:          : % e +/- % e" % (amp, sqrt(2./nt)*rms)
        print "Phase (ph) : % e +/- % e" % (ph, sqrt(2./nt)*rms/amp/(2.*pi))
        print "Phase (T0) : % e +/- % e" % (T0, sqrt(2./nt)*rms/amp/(2.*pi)/fbest)
        print "Offset     : % e +/- % e" % (offset, sqrt(1./nt)*rms)
        print 60*"-"
      else:
        print "Input points: %6d, frequency points: %6d" % (nt, len(self.freq))
        print 
        print "Maximum power   : %f " % (self.power[bbin])
        print blue("Best sine period") + ": %f +/- %f" % (1./fbest, Psin_err)
        # print 60*"-"

      
#    def do_stats(self):
#        """
#        Some statistics about the periodogram
#        """
#        # Index with maximum power
#        max_power_idx = argmax(self.power)
#        # Maximum (unnormalized) power
#        max_power = self._upower[max_power_idx]
#    
#        rms = sqrt(self._YY * (1. - max_power))
#        
#        # Get the curvature in the power peak by fitting a parabola y=aa*x^2
#        if (max_power_idx > 1) and (max_power_idx < len(self.freq)-2):
#            # Shift the parabola origin to the power peak
#            xh = (self.freq[max_power_idx-1:max_power_idx+2] - self.freq[max_power_idx])**2
#            yh = self._upow[max_power_idx-1:max_power_idx+2] - self._upow[max_power_idx]
#            # Calculate the curvature (final equation from least squares)
#            aa = sum(yh*xh)/sum(xh*xh)
#            nt = float(self.N)
#            f_err = sqrt(-2./nt * max_power/aa*(1.-max_power)/max_power)
#            Psin_err = sqrt(-2./nt* max_power/aa*(1.-max_power)/max_power) / self.freq[max_power_idx]**2
#        else:
#            f_err = None
#            Psin_err = None
#            raise ValueError("WARNING: Highest peak is at the edge of the frequency range.\nNo output of frequency error.\nIncrease frequency range to sample the peak maximum.")
#
#            
#        fbest = self.freq[max_power_idx]
#        amp = sqrt(self._a[max_power_idx]**2 + self._b[max_power_idx]**2)
#        ph  = arctan2(self._a[max_power_idx], self._b[max_power_idx]) / (2.*pi)
#        T0  = min(self.th) - ph/fbest
#        # Re-add the mean
#        offset = self._off[max_power_idx] + self._Y
#    
#        # Statistics
#        print self.name + " - statistical output"
#        print "-----------------------------------"
#        print "Number of input points:     %6d" % (nt)
#        print "Weighted mean of dataset:   % e" % (self._Y)
#        print "Weighted rms of dataset:    % e" % (sqrt(self._YY))
#        print "Time base:                  % e" % (max(self.th) - min(self.th))
#        print "Number of frequency points: %6d" % (len(self.freq))
#        print
#        print "Maximum power, p :    % e " % (self.power[max_power_idx])
#        print "Maximum power (without normalization):   %e" % (max_power)
#        print "Normalization    : ", self.norm
#        print "RMS of residuals :    % e " % (rms)
#        if self.error is not None:
#          print "  Mean weighted internal error:  % e" %(sqrt(nt/sum(1./self.error**2)))
#        print "Best sine frequency : % e +/- % e" % (fbest, f_err)
#        print "Best sine period    : % e +/- % e" % (1./fbest, Psin_err)
#        print "Amplitude:          : % e +/- % e" % (amp, sqrt(2./nt)*rms)
#        print "Phase (ph) : % e +/- % e" % (ph, sqrt(2./nt)*rms/amp/(2.*pi))
#        print "Phase (T0) : % e +/- % e" % (T0, sqrt(2./nt)*rms/amp/(2.*pi)/fbest)
#        print "Offset     : % e +/- % e" % (offset, sqrt(1./nt)*rms)
#        print "-----------------------------------"



class MCMC:
    def __init__(self, filename_glob, burnin=0):
        self.chain_filenames = sorted(glob.glob(filename_glob))
        self.nchains = len(self.chain_filenames)

        self.burnin = burnin

    def reset_burnin(self, burnin):
        self.burnin = burnin
        self.read_chains()

    def read_chains(self):
        nc = self.nchains
        # black magic to build input from file list while skipping headers
        finput = [FileInput(f) for f in sorted(self.chain_filenames)]
        iterables = [islice(f, self.burnin+1, None) for f in finput]
        files = chain(*iterables)

        chains = np.loadtxt(files, unpack=True)

        self.nplanets = (chains.shape[0]-3)/5
        self.nsamples = chains.shape[1] / nc  # number of samples per chain

        self.chains = np.array(chains)
        self.chains = np.ma.masked_invalid(self.chains) # remove NaNs

        parameter_names = ['period', 'K', 'ecc', 'omega', 't0']
        parameter_names = np.array([[par + '_'+str(i) for par in parameter_names] for i in range(self.nplanets)])

        p = namedtuple('parameter', parameter_names.flat)

        # this makes a namedtuple of the ith chain
        # self.trace = p._make(self.chains[i,2:-1,:])

        # this makes a namedtuple of all the chains together
        self.trace = p._make(self.chains[2:-1,:])

    def kde1d(self, parameter, ind=None, npoints=100):
        try:
            trace_dict = self.trace._asdict()
            t = trace_dict[parameter]
            gkde = scipy.stats.gaussian_kde(t)
            if ind is None:
                ind = np.linspace(min(t), max(t), npoints)
            return ind, gkde.evaluate(ind)
        except KeyError:
            print 'No parameter named', parameter
            print 'Try one of these:', self.trace._fields

    def kde2d(self, parameters, ind=None, npoints=100):


        """ Perform kernel density estimation on the joint distribution
        of two parameters.

        parameters: list with two parameters
        ind, optional: mesh grid on which to evaluate the kde
        """
        p = np.atleast_1d(parameters)
        assert np.shape(p)[0] == 2, 'Give me two parameters!'

        trace_dict = self.trace._asdict()
        for par in p:
            try: 
                t = trace_dict[par]
            except KeyError:
                print 'No parameter named', par
                print 'Try one of these:', self.trace._fields
                return

        print 'Estimating kde for', p[0], 'and', p[1]
        t1 = trace_dict[p[0]]
        t2 = trace_dict[p[1]]
        gkde = scipy.stats.gaussian_kde([t1, t2])
        x,y = np.mgrid[min(t1) : max(t1) : (max(t1)-min(t1)) / 20, 
                       min(t2) : max(t2) : (max(t2)-min(t2)) / 20]
        z = np.array(gkde.evaluate([x.flatten(),y.flatten()])).reshape(x.shape)
        return x, y, z

    def do_plot_some_solutions(self, system):
        """ Plot the observed radial velocities as a function of time together
        with a random choice of solutions from the MCMC
        Data from each file are color coded and labeled.
        """
        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = system.time, system.vrad, system.error # temporaries
        tt = system.get_time_to_plot()
        vel = np.zeros_like(tt)

        n_solutions = 20
        # random n_solutions indices from chains array
        ind = np.random.choice(self.chains[0,:].size, size=n_solutions, replace=False)
        # parameters corresponding to those random indices
        par = self.chains[2:, ind]

        ## best solution found
        best_ind = np.argmax(self.chains[1, :])
        par_best = self.chains[2:, best_ind]

        newFig=True        
        if newFig: 
            plt.figure()

        # plot each solution
        for i in range(n_solutions):
            p = par[:,i]
            args = [tt] + list(p) + [vel]
            get_rvn(*args)
            if i == 0:
                plt.plot(tt, vel, '-k', alpha=0.5, label='random')
            else:
                plt.plot(tt, vel, '-k', alpha=0.5)

        # plot best solution
        args = [tt] + list(par_best) + [vel]
        get_rvn(*args)
        plt.plot(tt, vel, '-g', lw=2.5, label='best')

        # plot each files' values
        for i, (fname, [n, nout]) in enumerate(sorted(system.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
            #                     height=err[:m], beam=0.5,\
            #                     pen=pg.mkPen(None))
                                # pen={'color': 0.8, 'width': 2})
            # p.addItem(e)
            # p.plot(t[:m], rv[:m], symbol='o')
            plt.errorbar(t[:m], rv[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=os.path.basename(fname))
            t, rv, err = t[m:], rv[m:], err[m:]
        
        plt.xlabel('Time [days]')
        plt.ylabel('RV [km/s]')
        plt.legend()
        plt.tight_layout()
        plt.ticklabel_format(useOffset=False)

    def do_plot_trace(self, parameter):

        if parameter not in self.trace._fields:
            print 'ERROR'
            return

        i = self.trace._fields.index(parameter)
        s = self.nsamples

        t = self.trace[i]

        newFig=True        
        if newFig: 
            plt.figure()

        for j in range(self.nchains):
            plt.plot(t[j*s : (j+1)*s])


    def do_triangle_plot(self):

        triangle_plot(self.chains[2:].T, quantiles=[0.5])






