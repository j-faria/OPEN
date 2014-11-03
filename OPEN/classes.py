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
from string import ascii_lowercase
from copy import copy
# other imports
import numpy as np
import matplotlib.pylab as plt
import scipy.stats
from scipy.signal import argrelmax
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
from .utils import day2year, rms, ask_yes_no, triangle_plot, triangle_plot_kde
from ext.get_rvN import get_rvn as get_rvn_os
from ext.get_rvN_MultiSite import get_rvn as get_rvn_ms
from ext.gp import gp_predictor
from ext.julian import caldat


# if the interpreter calls open.py from another directory, matplotlib will not
# access the matplotlibrc file and will revert to its defaults. This makes sure
# we use it!
# p = os.path.dirname(os.path.realpath(__file__)) # path to this file (classes.py)
# rc( rc_params_from_file( p+'/../matplotlibrc' ) )
# plt.ion()

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

        # by default
        self.units = 'km/s'

        # save the extra quantities as a namedtuple if we read them
        extra = namedtuple('Extra', data.keys(), verbose=False)
        self.extras = extra(**data)
        self.extras_names = data.keys()

        # time, vrad and error can change, 
        # the *_full ones correspond always to the full set
        self.time_full = copy(self.time)
        self.vrad_full = copy(self.vrad)
        self.error_full = copy(self.error)
        # same for extras
        self.extras_full = copy(self.extras)

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
        if len(self.provenance) > 1:
            stats += info + "Total number of measurements: %d\n" % sum([i[0] for i in self.provenance.values()])
            stats += sinfo + "Timespan : %f days = %f years   ---   %f JD, %f JD\n" % (tspan, day2year(tspan), min(t), max(t))
        else:
            stats += info + "Timespan : %f days = %f years   ---   %f JD, %f JD\n" % (tspan, day2year(tspan), min(t), max(t))
        stats += sinfo + "RV span  : %f km/s = %f m/s\n" % (rvspan, rvspan*1e3)
        stats += sinfo + "RV rms [m/s] : %f\n\n" % rms(rv*1e3)
        stats += sinfo + "{:14s} : {:10.3f}\n".format('<RV> [m/s]', np.mean(rv*1e3))
        if len(self.extras) == 0: 
            clogger.info(stats)
            return

        try:
            stats += sinfo + "{:14s} : {:10.3f}\n".format('<fwhm> [km/s]', np.mean(self.extras.fwhm))
        except AttributeError:
            pass
        try:
            stats += sinfo + "{:14s} : {:10.3f}\n".format('<contrast>', np.mean(self.extras.contrast) )
        except AttributeError:
            pass
        try:
            stats += sinfo + "{:14s} : {:10.3f}\n".format('<BIS> [km/s]', np.mean(self.extras.bis_span) )
        except AttributeError:
            pass
        try:            
            stats += sinfo + "{:14s} : {:10.3f}\n".format('<S_index> [MW]', np.mean(self.extras.s_mw) )
        except AttributeError:
            pass
        try:
            stats += sinfo + "{:14s} : {:10.3f}\n".format('<log(rhk)>', np.mean(self.extras.rhk))
        except AttributeError:
            pass
        finally:
            clogger.info(stats)

            

    def do_plot_obs(self, newFig=True, leg=True, save=None):
        """ Plot the observed radial velocities as a function of time.
        Data from each file are color coded and labeled.
        """

        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = self.time, self.vrad, self.error # temporaries
        
        if newFig: 
            fig = plt.figure()
        ax2 = fig.add_subplot(111)
        ax2.set_xlabel('Time [days]', labelpad=20)
        ax2.set_ylabel('RV [%s]'%self.units)
        
        ny, years = self.get_years_observations()
        ax1 = ax2.twiny()
        ax1.xaxis.tick_bottom()
        ax2.xaxis.tick_top()
        ax2.plot(years, np.ones_like(years), alpha=0) # Create a dummy plot

        # plot each files' values
        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
            #                     height=err[:m], beam=0.5,\
            #                     pen=pg.mkPen(None))
                                # pen={'color': 0.8, 'width': 2})
            # p.addItem(e)
            # p.plot(t[:m], rv[:m], symbol='o')
            ax1.errorbar(t[:m], rv[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err = t[m:], rv[m:], err[m:]
        
        if leg: ax1.legend()
        plt.tight_layout()
        ax2.ticklabel_format(useOffset=False)

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)

        # plt.show()
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
        
        poly = np.poly1d(drift)  # the polynomial evaluator
        p = np.polyval(drift, self.time)  # normal polynomial, for 2nd plot
        st = np.sort(self.time) # need this otherwise the plot gets messed up
        sp = np.polyval(drift, st) # "sorted" polynomial
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
        plt.ticklabel_format(useOffset=False)
        plt.show()

    def do_plot_extras(self, extra, save=None):
        """ Plot other observed quantities as a function of time.

        Parameters
        ----------
        extra: string
          One of the quantities available in system.extras
        save: None or string
          If not None, save figure with filename save (can have extension; default:png)
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
        elif extra == 'fwhm' and ask_yes_no('Should I remove a linear trend first? (y/N) ', False):
            m, b = np.polyfit(t, self.extras[i], 1)
            yp = np.polyval([m, b], t)
            plt.plot(t, self.extras[i]-yp, 'o', label=extra+' (linear trend removed)')
        else:
            plt.plot(t, self.extras[i], 'o', label=extra)
        
        plt.xlabel('Time [days]')
        plt.ylabel(extra + ' []')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05))
        plt.minorticks_on()
        plt.tight_layout()
        plt.show()

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)
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

        get_rvn = get_rvn_os
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
        # plt.figure()
        # gs = gridspec.GridSpec(keplerians, 1)
        # for i in range(keplerians):
        #   get_rvn(tt, P[i], K[i], ecc[i], omega[i], T0[i], gam, final)
        #   ax = plt.subplot(gs[i])
        #   # plot each file's values
        #   for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
        #       m = n-nout # how many values are there after restriction
        #       ax.errorbar(np.modf(abs(t[:m]-T0[i])/P[i])[0], rv[:m], yerr=err[:m], \
        #                    fmt='o'+colors[i], label=fname)
        #       t, rv, err = t[m:], rv[m:], err[m:]
        #   ax.plot( np.modf(abs(tt-T0[i])/P[i])[0], final, 'k-')

        plt.show()

    def do_plot_resid(self, newFig=True, leg=True, save=None):
        """ Plot the residuals from a fit as a function of time.
        Data from each file are color coded and labeled.
        """
        if self.fit is None:  # there is no fit yet
            clogger.warning(yellow('Warning: ')+'This system does not have a fit')
            return

        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, err = self.time, self.error # temporaries
        r = self.fit['residuals']
        
        if newFig: 
            fig = plt.figure()
        ax2 = fig.add_subplot(111)
        ax2.set_xlabel('Time [days]', labelpad=20)
        ax2.set_ylabel('Residuals [%s]'%self.units)
        
        ny, years = self.get_years_observations()
        ax1 = ax2.twiny()
        ax1.xaxis.tick_bottom()
        ax2.xaxis.tick_top()
        ax2.plot(years, np.ones_like(years), alpha=0) # Create a dummy plot

        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
            #                     height=err[:m], beam=0.5,\
            #                     pen=pg.mkPen(None))
                                # pen={'color': 0.8, 'width': 2})
            # p.addItem(e)
            # p.plot(t[:m], rv[:m], symbol='o')
            ax1.errorbar(t[:m], r[:m], yerr=err[:m], \
                         fmt='o'+colors[i], label=os.path.basename(fname))
            t, err = t[m:], err[m:]
        
        ax1.axhline(y=0, ls='--', color='k')

        if leg: ax1.legend()
        plt.tight_layout()
        ax2.ticklabel_format(useOffset=False)

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)
        # plt.show()
        # pg.QtGui.QApplication.exec_()

    def save_fit(self, params, chi2):
        """ Helper function to save results from a fit to the system """
        if self.fit is not None: # there is already a fit
            clogger.warning(yellow('Warning: ')+'Replacing older fit')

        self.fit = {}
        self.fit['params'] = params

        P, K, ecc, omega, T0, gam = [params[i::6] for i in range(6)]
        gam = gam[0]

        final = np.zeros_like(self.time)
        get_rvn = get_rvn_os
        get_rvn(self.time, P, K, ecc, omega, T0, gam, final)
        self.fit['residuals'] = self.vrad - final
        self.fit['chi2'] = sum((self.fit['residuals'] / self.error)**2)
        self.fit['chi2_reduced'] = self.fit['chi2'] / (len(self.time) - len(self.fit['params']))

        self.do_fit_stats()

    def do_fit_stats(self):
        """ Calculate and print some statistics related to the (current) fit """
        fit = self.fit

        nobs = len(self.time)
        try:
            k_trend = len(self.model['drift'])
        except (TypeError, KeyError):  # it raises a TypeError because self.model is None
            k_trend = 0
        k = len(fit['params']) + k_trend  # number of free parameters
        r = fit['residuals'] / self.error  # error weighted residuals

        # log-likelihood
        ll = 0.5 * nobs * np.log(1./2.*np.pi) - 0.5 * np.sum(r**2)
        aic = -2*ll + 2.*k
        aicc = aic + (2.*k*(k+1)) / (nobs-k-1.)
        bic = -2*ll + (k * np.log(nobs))

        self.fit['ll'] = ll
        self.fit['aic'] = aic
        self.fit['aicc'] = aicc
        self.fit['bic'] = bic

    # def do_stats(fit, x, y):
    #     def ll_aic_bic(fit):
    #         """
    #         Calculate model log-likelihood, AIC, and BIC information criteria
    #         """
    #         # ll = -(nobs*1/2)*(1+np.log(2*np.pi)) - (nobs/2.)*np.log(np.dot(e,e)/nobs)
    #         # print ll
    #         ll = 0.5 * nobs * np.log(1./2*np.pi) - 0.5 * np.sum(e**2)
    #         k = fit.nvarys
    #         aic = -2*ll + 2.*k
    #         aicc = aic + (2.*k*(k+1)) / (nobs-k-1.)
    #         bic = -2*ll + (fit.nvarys * np.log(nobs))
    #         # print ll, aic, bic
    #         return ll, aic, bic, aicc

    #     e = fit.residual
    #     nobs = fit.ndata

    #     # durbin watson statistic -- http://en.wikipedia.org/wiki/Durbin%E2%80%93Watson_statistic
    #     fit.DurbinWatson = sum(ediff1d(fit.residual)**2) / sum(fit.residual**2)

    #     # r squared, coefficient of determination
    #     fit.R2 = 1 - sum(fit.residual**2) / sum((y - mean(y))**2)
    #     fit.R2adj = 1 - (1-fit.R2)* ((nobs-1)/(nobs-fit.nvarys-1))   # adjusted R-square

    #     # log-likelihood, AIC, BIC
    #     fit.ll, fit.aic, fit.bic, fit.aicc = ll_aic_bic(fit)

    #     # residual skewness, kurtosis, and JB test for normality
    #     fit.skew = sp.stats.skew(e)
    #     fit.kurtosis = 3. + sp.stats.kurtosis(e)
    #     fit.JB = (e.shape[0] / 6.) * (fit.skew**2 + (1 / 4.) * (fit.kurtosis-3)**2)
    #     fit.JBpv = sp.stats.chi2.sf(fit.JB,2)


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

    def get_years_observations(self):
        if self.time[0] / 24e5 > 1:
            bjd0 = 0
        else:
            bjd0 = 24e5
        m1, _, y1 = caldat(bjd0+self.time[0])
        m2, _, y2 = caldat(bjd0+self.time[-1])
        # print m1, y1
        # print m2, y2
        # return len(np.arange(y1, y2)), np.arange(y1, y2)
        r = np.arange(y1,y2+1)
        # print min(r), max(r)
        return len(r), r

    def get_time_to_plot(self):
        """
        Returns sufficiently resolved time vector for plots
        """
        std = self.time.std()
        N = len(self.time)
        minim = self.time.min() - 1.*std
        maxim = self.time.max() + 1.*std
        return np.linspace(minim, maxim, 10*N)

    def is_in_extras(self, extra):
        return extra in self.extras._fields


class BasicTimeSeries():
    time = None
    vrad = None
    error = None


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

    def _plot(self, doFAP=False, faps=None, verts=None, newFig=True, axes=None, save=None, **kwargs):
        """
        Plot this periodogram.
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

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)

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


from time import strftime
from zipfile import ZipFile     

class MCMC_dream:
    """
    Base class for MCMC analysis, adapted for DREAM output
    """
    def __init__(self, filename_glob, burnin=0):
        if filename_glob.endswith('.zip'):
            with ZipFile(filename_glob) as zf:
                self.chain_filenames = [n for n in zf.namelist() if 'chain' in n]
                zf.extractall()
            
        else:
            self.chain_filenames = sorted(glob.glob(filename_glob))

        self.nchains = len(self.chain_filenames)
        self.burnin = burnin

        self.get_rvn = get_rvn_os

        self.read_chains()

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

    def compress_chains(self):
        tstr = strftime("%Y%m%d-%H%M%S")
        zfilename = 'demc-out-'+tstr+'.zip'

        with ZipFile(zfilename, 'a', compression=8) as zf:
            for f in self.chain_filenames:
                zf.write(f)
            zf.write('OPEN/demc/namelist1')

        print 'done!'


    def print_best_solution(self):
        ## output best solution
        best_ind = np.argmax(self.chains[1, :])
        par_best = self.chains[2:, best_ind]
        print # newline
        msg = yellow('RESULT: ') + 'Best solution is'
        clogger.info(msg)
        ## loop over planets
        print("%3s %12s %10s %10s %10s %15s %9s" % \
            ('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
        for i, planet in enumerate(list(ascii_lowercase)[:self.nplanets]):
            P, K, ecc, omega, T0, gam = [par_best[j::6] for j in range(6)]
            print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % 
                                (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )


    def confidence_intervals(self):
        for i, name in enumerate(self.trace._fields):
            print '%10s' % name,
            print '  ', 'mean:', np.mean(self.trace[i]), 'std:', np.std(self.trace[i])


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
            self.get_rvn(*args)
            if i == 0:
                plt.plot(tt, vel, '-k', alpha=0.5, label='random')
            else:
                plt.plot(tt, vel, '-k', alpha=0.5)

        # plot best solution
        args = [tt] + list(par_best) + [vel]
        self.get_rvn(*args)
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
        plt.show()

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

    def do_plot_best(self, system):
        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = system.time, system.vrad, system.error # temporaries
        tt = system.get_time_to_plot()
        vel = np.zeros_like(tt)

        ## best solution found
        best_ind = np.argmax(self.chains[1, :])
        par_best = self.chains[2:, best_ind]

        newFig=True        
        if newFig: 
            plt.figure()

        # plot best solution
        args = [tt] + list(par_best) + [vel]
        self.get_rvn(*args)
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
        plt.show()

    def do_triangle_plot(self):

        triangle_plot(self.chains[2:].T, quantiles=[0.5])


    def do_triangle_plot_kde(self):

        triangle_plot_kde(self.chains[2:].T)



class MCMC_nest:
    """
    Base class for MCMC analysis, adapted for MultiNest output
    """
    def __init__(self, root, context='111'):
        self.root = root
        with open(root+'stats.dat') as f:
            for line in f: pass  # loop over the file
            last = line  # to get the last line

        self.context = context

        # this seems like the easiest way to get the number of parameters
        self.npar = int(last.split()[0])

        self.gp = self.gp_only = False
        if (context[0] in ('2', '3')):
            self.gp = True
            self.read_gp_file()
            if context[0] == '3': 
                self.gp_only = True

        self.read_stats_file()

    def read_stats_file(self):
        """ 
        Read stats.dat file for mean, ML and MAP parameters
        """
        with open(self.root+'stats.dat') as f:
            stats = f.readlines()

        npar = self.npar
        self.nplanets = npar/5
        if self.gp and (npar == 5): 
            self.nplanets = 0
            self.gp_only = True

        try:
            self.NS_lnE = float(stats[0].split()[-3])
            self.NS_lnE_error = float(stats[0].split()[-1])
        except ValueError:
            self.NS_lnE, self.NS_lnE_error = 0, 0
        try:
            self.INS_lnE = float(stats[1].split()[-3])
            self.INS_lnE_error = float(stats[1].split()[-1])
        except ValueError:
            self.INS_lnE, self.INS_lnE_error = 0, 0

        ## mean params
        self.par_mean = [float(s.split()[1]) for s in stats[4:4+npar]]
        # P_mean, K_mean, ecc_mean, omega_mean, t0_mean, vsys_mean = par_mean
        self.par_sigma = [float(s.split()[2]) for s in stats[4:4+npar]]
        # P_sigma, K_sigma, ecc_sigma, omega_sigma, t0_sigma, vsys_sigma = par_sigma

        ## MLE
        start, end = 4+npar+3, 4+npar+3+npar
        self.par_mle = [float(s.split()[1]) for s in stats[start:end]]
        # P_mle, K_mle, ecc_mle, omega_mle, t0_mle, vsys_mle = par_mle

        ## MAP
        start, end = 4+2*3+2*npar, 4+2*3+3*npar
        self.par_map = [float(s.split()[1]) for s in stats[start:end]]
        # P_map, K_map, ecc_map, omega_map, t0_map, vsys_map = par_map


    def read_gp_file(self):
        filename = self.root+'gp.dat'
        try:
            self.pred_t, self.pred_y, self.pred_std = np.loadtxt(filename, unpack=True)
        except IOError:
            self.pred_t, self.pred_y, self.pred_std = (np.nan, np.nan, np.nan)

    def read_iterations(self):
        """
        Read the [root].txt file.
        Columns have sample probability, -2*loglikehood, parameter values. 
        Sample probability is the sample prior mass multiplied by its 
        likelihood and normalized by the evidence.
        """
        filename = self.root + '.txt'

        data = np.loadtxt(filename, unpack=True)
        self.posterior_samples = data

        parameter_names = ['period', 'K', 'ecc', 'omega', 't0']
        parameter_names = np.array([[par + '_'+str(i) for par in parameter_names] for i in range(self.nplanets)])

        p = namedtuple('parameter', parameter_names.flat)

        # this makes a namedtuple of all the parameter values
        self.trace = p._make(data[2:-1,:])

    def compress_chains(self):
        tstr = strftime("%Y%m%d-%H%M%S")
        zfilename = 'demc-out-'+tstr+'.zip'

        with ZipFile(zfilename, 'a', compression=8) as zf:
            for f in self.chain_filenames:
                zf.write(f)
            zf.write('OPEN/demc/namelist1')

        print 'done!'


    def print_best_solution(self):
        ## output best solution
        best_ind = np.argmax(self.chains[1, :])
        par_best = self.chains[2:, best_ind]
        print # newline
        msg = yellow('RESULT: ') + 'Best solution is'
        clogger.info(msg)
        ## loop over planets
        print("%3s %12s %10s %10s %10s %15s %9s" % \
            ('', 'P[days]', 'K[km/s]', 'e', unichr(0x3c9).encode('utf-8')+'[deg]', 'T0[days]', 'gam') )
        for i, planet in enumerate(list(ascii_lowercase)[:self.nplanets]):
            P, K, ecc, omega, T0, gam = [par_best[j::6] for j in range(6)]
            print("%3s %12.1f %10.2f %10.2f %10.2f %15.2f %9.2f" % 
                                (planet, P[i], K[i], ecc[i], omega[i], T0[i], gam[i]) )


    def confidence_intervals(self):
        for i, name in enumerate(self.trace._fields):
            print '%10s' % name,
            print '  ', 'mean:', np.mean(self.trace[i]), 'std:', np.std(self.trace[i])


    def save_fit_to(self, system):
        
        if system.fit is not None: # there is already a fit
            clogger.warning(yellow('Warning: ')+'Replacing older fit')

        system.fit = {}
        
        t, rv, err = system.time, system.vrad, system.error # temporaries
        vel = np.zeros_like(t)

        nobserv = len(system.provenance)  # number of observatories
        if nobserv > 1:
            vsys = self.par_map[-nobserv+1:]
            get_rvn = get_rvn_ms
        else:
            vsys = self.par_map[-1]
            get_rvn = get_rvn_os

        if self.npar == 1:  # systematic velocity only
            vel = self.par_map[0]

        ## MAP estimate of the parameters
        elif self.gp:
            par_map = self.par_map[:-4] 
            hyper_map = self.par_map[-4:]
            if self.gp_only:
                vel = gp_predictor(t, rv, err, par_map, hyper_map, 'constant')
            else:
                vel = gp_predictor(t, rv, err, par_map, hyper_map, 'keplerian')
        else:
            par_map = self.par_map

            P = par_map[:-1:5]
            K = par_map[1:-1:5]
            ecc = par_map[2:-1:5]
            omega = par_map[3:-1:5]
            t0 = par_map[4:-1:5]
            par = [P, K, ecc, omega, t0, vsys]

            args = [t] + par + [vel]
            get_rvn(*args)


        system.fit['params'] = self.par_map
        system.fit['residuals'] = system.vrad - vel



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
        plt.ylabel('RV [%s]'%system.units)
        plt.legend()
        plt.tight_layout()
        plt.ticklabel_format(useOffset=False)
        plt.show()

    def do_plot_trace(self, parameter):

        if parameter not in self.trace._fields:
            print 'ERROR'
            return

        i = self.trace._fields.index(parameter)
        t = self.trace[i]

        newFig=True        
        if newFig: 
            plt.figure()

        plt.plot(t)
        plt.ylabel(self.trace._fields[i])

    def do_plot_map(self, system):
        colors = 'kbgrcmyw' # lets hope for less than 9 data-sets
        t, rv, err = system.time, system.vrad, system.error # temporaries
        tt = system.get_time_to_plot()
        vel = np.zeros_like(tt)

        nobserv = len(system.provenance)  # number of observatories
        if nobserv > 1:
            vsys = self.par_map[-nobserv+1:]
            get_rvn = get_rvn_ms
        else:
            vsys = self.par_map[-1]
            get_rvn = get_rvn_os

        ## MAP estimate of the parameters
        if self.npar == 1:  # systematic velocity only
            vel = self.par_map[0] * np.ones_like(tt)
        elif self.gp:
            par_map = self.par_map[:-4] 
            hyper_map = self.par_map[-4:]
            if self.gp_only:
                pred = gp_predictor(t, rv, err, par_map, hyper_map, 'constant')
            else:
                pred = gp_predictor(t, rv, err, par_map, hyper_map, 'keplerian')
        else:
            par_map = self.par_map

        newFig=True
        if newFig:
            plt.figure()

        gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1, sharey=ax1)

        if self.npar == 1:  # systematic velocity only
            ax1.plot(tt, vel, '-g', lw=2.5, label='MAP')
        elif not self.gp_only:
            # plot best solution
            P = par_map[:-1:5]
            K = par_map[1:-1:5]
            ecc = par_map[2:-1:5]
            omega = par_map[3:-1:5]
            t0 = par_map[4:-1:5]
            par = [P, K, ecc, omega, t0, vsys]
            # print par

            if nobserv > 1:
                observ = np.ones_like(vel)
                args = [tt] + par + [observ] + [vel]
            else:
                args = [tt] + par + [vel]
            get_rvn(*args)
            ax1.plot(tt, vel, '-g', lw=2.5, label='MAP')

        # plot GP predictions
        if self.gp:
            ax1.plot(self.pred_t, self.pred_y, '-k', lw=1.5, label='GP mean')
            ax1.fill_between(self.pred_t, y1=self.pred_y-2*self.pred_std, y2=self.pred_y+2*self.pred_std,
                                          color='k', alpha=0.3, label='2*std')

        vel = np.zeros_like(t)
        if self.gp_only or self.npar == 1:
            vel = self.par_map[0] * np.ones_like(tt)
        else:
            args = [t] + par + [vel]
            get_rvn(*args)

        for i, (fname, [n, nout]) in enumerate(sorted(system.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
            #                     height=err[:m], beam=0.5,\
            #                     pen=pg.mkPen(None))
                                # pen={'color': 0.8, 'width': 2})
            # p.addItem(e)
            # p.plot(t[:m], rv[:m], symbol='o')

            # plot each files' values
            ax1.errorbar(t[:m], rv[:m], yerr=err[:m], fmt='o'+colors[i], label=os.path.basename(fname))
            # plot residuals
            if self.gp:
                ax2.errorbar(t[:m], rv[:m]-pred[:m], yerr=err[:m], fmt='o'+colors[i], label=fname)
            else:
                ax2.errorbar(t[:m], rv[:m]-vel[:m], yerr=err[:m], fmt='o'+colors[i], label=fname)

            t, rv, err = t[m:], rv[m:], err[m:]


        ax2.axhline(y=0, ls='--', color='k')

        ax2.set_xlabel('Time [days]')
        ax1.set_ylabel('RV [%s]'%system.units)
        ax2.set_ylabel('Residuals [%s]'%system.units)
        ax1.legend()
        plt.tight_layout()
        plt.ticklabel_format(useOffset=False)
        plt.show()


    def do_plot_map_phased(self, system, noname=False, plot_gp=True):
        get_rvn = get_rvn_os
        colors = 'kbgrcmyw' # lets hope for less than 9 data-sets
        t, rv, err = system.time, system.vrad, system.error # temporaries
        tt = system.get_time_to_plot()
        vel = np.zeros_like(tt)

        if self.gp:
            self.npar

        ## MAP estimate of the parameters
        if self.gp:
            par_map = self.par_map[:-4]  # don't care about the hyperparameters for now
        else:
            par_map = self.par_map
            if self.npar in (7, 12):
                par_map.pop(-2)

        newFig=True
        if newFig:
            plt.figure()

        # plot best solution
        P = par_map[:-1:5]
        K = par_map[1:-1:5]
        ecc = par_map[2:-1:5]
        omega = par_map[3:-1:5]
        t0 = par_map[4:-1:5]
        par = [P, K, ecc, omega, t0, par_map[-1]]

        args = [tt] + par + [vel]
        get_rvn(*args)
        phase = ((tt - t0) / P) % 1.0
        plt.plot(np.sort(phase), vel[np.argsort(phase)], '-g', lw=2.5, label='MAP')
        plt.plot(np.sort(phase)+1, vel[np.argsort(phase)], '-g', lw=2.5)
        plt.plot(np.sort(phase)-1, vel[np.argsort(phase)], '-g', lw=2.5)

        # plot GP predictions
        if self.gp and plot_gp:
            phase = ((self.pred_t - t0) / P) % 1.0
            indices = np.argsort(phase)
            plt.plot(np.sort(phase), self.pred_y[indices], '-k', lw=0.5, alpha=0.6, label='GP mean')
            # plt.fill_between(np.sort(phase), 
            #                  y1=self.pred_y[indices]-2*self.pred_std[indices], 
            #                  y2=self.pred_y[indices]+2*self.pred_std[indices],
            #                  color='k', alpha=0.3, label='2*std')

        # plot each files' values
        for i, (fname, [n, nout]) in enumerate(sorted(system.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
            #                     height=err[:m], beam=0.5,\
            #                     pen=pg.mkPen(None))
                                # pen={'color': 0.8, 'width': 2})
            # p.addItem(e)
            # p.plot(t[:m], rv[:m], symbol='o')
            phase = ((t[:m] - t0) / P) % 1.0
            plt.errorbar(np.sort(phase), rv[np.argsort(phase)], yerr=err[np.argsort(phase)], \
                         fmt='o'+colors[i], label=os.path.basename(fname))
            plt.errorbar(np.sort(phase)+1, rv[np.argsort(phase)], yerr=err[np.argsort(phase)], \
                         fmt='o'+colors[i])
            plt.errorbar(np.sort(phase)-1, rv[np.argsort(phase)], yerr=err[np.argsort(phase)], \
                         fmt='o'+colors[i])
            t, rv, err = t[m:], rv[m:], err[m:]
        
        plt.xlim([-0.2, 1.2])
        plt.xlabel('Phase')
        plt.ylabel('RV [%s]'%system.units)
        if not noname: plt.legend()
        plt.tight_layout()
        plt.ticklabel_format(useOffset=False)
        plt.show()


    def do_triangle_plot(self):

        triangle_plot(self.posterior_samples[2:,:].T, quantiles=[0.5])


    def do_triangle_plot_kde(self):

        triangle_plot_kde(self.posterior_samples[2:,:].T)



