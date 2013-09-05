# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
import numpy
import rvIO
from utils import unique
from logger import clogger, logging

import matplotlib.pylab as plt
import plots_config


class rvSeries:
    """
    A container class that holds the observed radial velocity curve. 
    
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

    def __init__(self, *filenames):
        
        assert len(filenames)>=1, "Need at least one file to read"
        # don't repeat files
        filenames = unique(filenames)

        t, rv, err, self.provenance = rvIO.read_rv(*filenames, verbose=False)
        self.time, self.vrad, self.error = (t, rv, err)
        # time, vrad and error can change, 
        # the *_full ones correspond always to the full set
        self.time_full = self.time
        self.vrad_full = self.vrad
        self.error_full = self.error

    # associated model to be adjusted to the data
    # this will be a dict with the following key:value pairs:
    #   k : number of keplerians in the model
    #   d : degree of drift in the model
    #   drift : array with parameters from drift fit (the polynomial coefficients)
    model = None

    def do_plot_obs(self):
        """ Plot the observed radial velocities as a function of time.
        Data from each file is color coded and labeled.
        """

        colors = 'rgbmk' # possible colors
        t, rv, err = self.time, self.vrad, self.error # temporaries
        
        plt.figure()
        # plot each files' values
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            print fname, m

            plt.errorbar(t[:m], rv[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err = t[m:], rv[m:], err[m:]
        
        plt.xlabel('Time [days]')
        plt.ylabel('RV [m/s]')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def do_plot_drift(self):
        """ Plot the observed radial velocities as a function of time, plus an
        extra drift of specified degree (see *mod*). Lower panel presents RV 
        minus drift. Data from each file is color coded and labeled.
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
        ax2 = plt.subplot(2,1,2, sharex=ax1, sharey=ax1)
        # plot each file's values minus the drift
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            ax2.errorbar(t[:m], rv[:m]-p[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err, p = t[m:], rv[m:], err[m:], p[m:]

        plt.xlabel('Time [days]')
        plt.ylabel('RV [m/s]')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def get_nyquist(self):
        """
        Calculate the average Nyquist frequency.
          
        Returns
        -------
        Nyquist frequency : float
            Half the sampling frequency of the time series.
        """
        return 0.5*1./numpy.mean(self.time[1::]-self.time[0:-1])
        
        
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

