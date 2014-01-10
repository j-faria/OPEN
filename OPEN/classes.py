# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

# other imports
import numpy
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from collections import namedtuple # this requires Python >= 2.6

# intra-package imports
import rvIO
from .utils import unique
from .logger import clogger, logging
from ext.get_rvN import get_rvn


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

        assert len(filenames)>=1, "Need at least one file to read"
        # don't repeat files
        filenames = unique(filenames)

        # verbose?
        try:
          verbose = kwargs.pop('verbose')
        except KeyError:
          verbose = False
          
        # skip header lines? (no by default)
        try:
          skip = kwargs.pop('skip')
        except KeyError:
          skip = 0

        # in which format are the files? (drs35 by default)
        try:
          format = kwargs.pop('format').lower()
        except (KeyError, AttributeError): # this way we also catch the default None
          format = 'drs35'

        # read data
        try:
          t, rv, err, self.provenance, extras = \
               rvIO.read_rv(*filenames, verbose=verbose, skip=skip, format=format)
        except ValueError:
          from shell_colors import red
          msg = red('ERROR: ') + 'If your files have header lines set --skip option.\n'
          clogger.fatal(msg)
          return

        self.time, self.vrad, self.error = (t, rv, err)

        # save the extra quantities as a namedtuple if we read them
        if format == 'drs35': # default
          extras_names = ['fwhm', 'contrast', 'bis_span', 'noise', 
                          's_mw', 'sig_s', 'rhk', 'sig_rhk', 'sn_CaII', 
                          'sn10', 'sn50', 'sn60']
        elif format == 'drs34' or format == 'coralie':
          extras_names = ['fwhm', 'contrast', 'bis_span', 'noise', 'sn10', 'sn50', 'sn60']
        else:
          extras_names = []
        extra = namedtuple('Extra', extras_names, verbose=False)
        self.extras = extra._make(extras)
        self.extras_names = extras_names

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
    #   chi2 : reduced(!) chi square value of the fit
    fit = None

    def do_plot_obs(self):
        """ Plot the observed radial velocities as a function of time.
        Data from each file are color coded and labeled.
        """
        # import pyqtgraph as pg

        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = self.time, self.vrad, self.error # temporaries
        
        plt.figure()
        # p = pg.plot()
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
        plt.legend()
        plt.tight_layout()
        plt.show()
        # pg.QtGui.QApplication.exec_()

    def do_plot_drift(self):
        """ Plot the observed radial velocities as a function of time, plus an
        extra drift of specified degree (see *mod*). Lower panel presents RV 
        minus drift. Data from each file are color coded and labeled.
        """
        colors = 'rgbmk' # possible colors
        t, rv, err = self.time, self.vrad, self.error # temporaries
        
        plt.close(1)
        plt.figure(1)

        ax1 = plt.subplot(2,1,1)
        # plot each file's values
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            ax1.errorbar(t[:m], rv[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err = t[m:], rv[m:], err[m:]
        
        drift = self.model['drift']
        p = numpy.polyval(drift, self.time) # normal polynomial, for 2nd plot
        st = numpy.sort(self.time) # need this otherwise the plot gets messed up
        sp = numpy.polyval(drift, st) # "sorted" polynomial
        ax1.plot(st, sp, 'y-')

        t, rv, err = self.time, self.vrad, self.error # temporaries
        ax2 = plt.subplot(2,1,2, sharex=ax1)
        # plot each file's values minus the drift
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            ax2.errorbar(t[:m], rv[:m]-p[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err, p = t[m:], rv[m:], err[m:], p[m:]

        plt.xlabel('Time [days]')
        plt.ylabel('RV [km/s]')
        # plt.legend()
        plt.tight_layout()
        plt.show()

    def do_plot_extras(self, extra):
        """ Plot other observed quantities as a function of time.

        Parameters
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

        tt = numpy.linspace(self.time.min()-200, self.time.max()+300, 400)
        final = numpy.zeros_like(tt)
        
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

        # redo this...
        t, rv, err = self.time, self.vrad, self.error # temporaries
        final = numpy.zeros_like(t)
        get_rvn(t, P, K, ecc, omega, T0, gam, final)

        ax2 = plt.subplot(gs[1], sharex=ax1)
        # plot residuals
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            ax2.errorbar(t[:m], rv[:m]-final[:m], yerr=err[:m], \
                       fmt='o'+colors[i], label=fname)
            t, rv, err, final = t[m:], rv[m:], err[m:], final[m:]
        ax2.axhline(y=0, xmin=0, xmax=1, ls='--', color='k')

        # # redo this...
        # t, rv, err = self.time-self.time.min(), self.vrad, self.error # temporaries
        # final = numpy.zeros_like(tt)

        # ### phased RV curves
        # plt.figure()
        # gs = gridspec.GridSpec(keplerians, 1)
        # for i in range(keplerians):
        #   get_rvn(tt, P[i], K[i], ecc[i], omega[i], T0[i], gam, final)
        #   ax = plt.subplot(gs[i])
        #   # plot each file's values
        #   for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
        #       m = n-nout # how many values are there after restriction
        #       ax.errorbar(t[:m]/P[i], rv[:m], yerr=err[:m], \
        #                    fmt='o'+colors[i], label=fname)
        #       t, rv, err = t[m:], rv[m:], err[m:]
        #   ax.plot((tt-min(tt))/P[i], final, 'k-')

        plt.show()

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
            return 0.5*1./numpy.min(self.time[1::]-self.time[0:-1])
        else:
            return 0.5*1./numpy.mean(self.time[1::]-self.time[0:-1])

    def get_time_to_plot(self):
        """
        Returns sufficiently resolved time vector for plots
        """
        std = self.time.std()
        N = len(self.time)
        minim = self.time.min() - 1.*std
        maxim = self.time.max() + 1.*std
        return numpy.linspace(minim, maxim, 2*N)
        
        
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
        temp = sorted(zip(self.power, self.freq))

        # largest peak
        fmax1 = (temp[nf-1])[1]
        if n==1:
            if output_period: 
                return fmax1, 1./fmax1
            else:
                return fmax1

        # second largest
        fmax2 = fmax1
        i = 1
        while abs(fmax2 - fmax1) < 0.01:
            i += 1
            fmax2 = (temp[nf-i])[1]

        if n==2:
            if output_period: 
                return [(fmax1, 1./fmax1), 
                        (fmax2, 1./fmax2)]
            else:
                return fmax1, fmax2

        # third largest
        fmax3 = fmax2
        j = i
        while abs(fmax3 - fmax2) < 0.01 or abs(fmax3 - fmax1) < 0.01:
            j += 1
            fmax3 = (temp[nf-j])[1]

        if n==3:
            if output_period: 
                return [(fmax1, 1./fmax1), 
                        (fmax2, 1./fmax2), 
                        (fmax3, 1./fmax3)]
            else:
                return fmax1, fmax2, fmax3

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

    def _plot(self, doFAP=None):
      """
        Create a plot.
      """
      xlabel = 'Period [d]'
      ylabel = 'Power'

      self.fig = plt.figure()
      self.ax = self.fig.add_subplot(1,1,1)
      self.ax.set_title("Normalized periodogram")
      self.ax.set_xlabel(xlabel)
      self.ax.set_ylabel(ylabel)
      self.ax.semilogx(1./self.freq, self.power, 'b-')
      # plot FAPs
      if doFAP is None: 
        pass
      elif doFAP is True: # do default FAPs of 10%, 1% and 0.1%
        pmin = 1./self.freq.min()
        pmax = 1./self.freq.max()
        plvl1 = self.powerLevel(0.1) # 10% FAP
        plvl2 = self.powerLevel(0.01) # 1% FAP
        plvl3 = self.powerLevel(0.001) # 0.1% FAP
        self.ax.semilogx([pmin, pmax],[plvl1, plvl1],'k-')
        self.ax.semilogx([pmin, pmax],[plvl2, plvl2],'k--')
        self.ax.semilogx([pmin, pmax],[plvl3, plvl3],'k:')

      plt.tight_layout()
      plt.show()
      # p = pg.plot(1./self.freq, self.power, title="Periodogram")
      # p.plotItem.setLogMode(1./self.freq, self.power)
      # pg.QtGui.QApplication.exec_()

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

