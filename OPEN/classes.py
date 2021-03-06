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
from .logger import clogger, logging
from .shell_colors import yellow, blue
from .utils import unique, get_tt, day2year, rms, ask_yes_no, triangle_plot, triangle_plot_kde, get_star_name
from ext.get_rvN import get_rvn
# from ext.get_rvN_MultiSite import get_rvn as get_rvn_ms
from ext.gp import gp_predictor
from ext.julian import caldat


# matplotlib parameters for publication plots
import matplotlib.ticker as ticker
import re
# plt.rcParams['text.latex.preamble'].append(r'\usepackage{lmodern}')
# plt.rcParams["text.latex.preamble"].append(r'\mathchardef\mhyphen="2D')
params = {'text.latex.preamble': [r'\usepackage{lmodern}', 
                                  r'\usepackage{amsfonts,amsmath,amssymb}',
                                  r'\mathchardef\mhyphen="2D',
                                  r'\DeclareMathOperator{\ms}{m\,s^{-1}}'
                                  r'\DeclareMathOperator{\kms}{km\,s^{-1}}'],
          'text.usetex' : True,
          'font.size'   : 8,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          'axes.unicode_minus': True,
          }

class MyFormatter(ticker.ScalarFormatter):
    """ Axis ticks formatter to replace the big minus sign with a smaller, prettier one. """
    def __call__(self, x, pos=None):
        # call the original LogFormatter
        rv = ticker.ScalarFormatter.__call__(self, x, pos)
        # check if we really use TeX
        if plt.rcParams["text.usetex"]:
            # if we have the string ^{- there is a negative exponent
            # where the minus sign is replaced by the short hyphen
            rv = re.sub(r'-', r'\mhyphen', rv)
        return rv


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
        skip = kwargs.get('skip', 0)

        # verbosity (on by default)
        verbose = kwargs.get('verbose', True)

        # read data
        try:
          data, self.provenance = rvIO.read_rv(*filenames, skip=skip, verbose=verbose)
        except ValueError as e:
          print e
          msg = red('ERROR: ') + 'If your files have header lines set --skip option.\n'
          clogger.fatal(msg)
          return

        # # the following is the most readable and less nested way I found of doing this
        # def get_time(d):
        #     try:
        #         return d.pop('jdb')
        #     except KeyError:
        #         pass
        #     try:
        #         return d.pop('time')
        #     except KeyError:
        #         pass
        #     raise IOError('The .rdb file header should have standard names: jdb/time, rv/vrad, err/svrad.')

        # def get_RV(d):
        #     try:
        #         return d.pop('vrad'), d.pop('svrad')
        #     except KeyError:
        #         pass
        #     try:
        #         return d.pop('rv'), d.pop('err')
        #     except KeyError:
        #         pass
        #     raise IOError('The .rdb file header should have standard names: jdb/time, rv/vrad, err/svrad.')

        # self.time = get_time(data)
        # self.vrad, self.error = get_RV(data)

        self.data = data
        
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

        # each time a fit is saved to the system, the residuals should be copied
        # to the following array. We assume less than 10 planets will be removed
        self.resid = np.zeros((len(self.time), 10))
        # the "first" residuals are the original radial velocities
        self.resid[:,0] = self.vrad_full

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

    # the number of planets (signals) whose residuals have been calculated
    # after removing a planet, update this and copy residuals to the array
    # self.resid[:, self.n_removed_planets]
    n_removed_planets = 0

    # properties of the star
    star_mass = 1.0  # in Msun, default 1.0


    @property
    def star_name(self):
        return self.__star_name
    @star_name.getter
    def star_name(self):
        self.__star_name = get_star_name(self)
        return self.__star_name

    def save(self, filename, **kwargs):
        rvIO.write_rv(self, filename, **kwargs)

    def stats(self):
        # verbose stats about data
        t, rv = self.time, self.vrad
        info = blue('INFO: ') 
        sinfo = blue('    : ') 
        stats = None
        tspan = max(t) - min(t)
        self.timespan = tspan  # in days!
        self.timespan_years = day2year(tspan)  # in years!
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

            
    def set_errors_indicators(self):
        """
        Set the photon noise uncertainty in various diagnoses,
        using the scaling coefficients (for HARPS) from Santerne et al. (2015, Table 1).
        It only changes (and adds to) the `extras` quantities.
        """
        if self.is_in_extras('noise'):  # we can only scale from the RV photon noise
            noise = self.extras.noise

            noiseFWHM = 2. * noise
            noiseBIS = 2. * noise

            self.data['sig_fwhm'] = noiseFWHM
            self.data['sig_bis'] = noiseBIS

            ### here we should check if other indicators are present
            ### BiGauss = 2.1
            ### Vspan = 1.5
            ### Wspan = 5.8
            ### Vasy = 5.6e10  (do not use according to Santerne+2015)

            extra = namedtuple('Extra', self.data.keys(), verbose=False)
            self.extras = extra(**self.data)
            self.extras_names = self.data.keys()

        else:
            from shell_colors import red
            msg = red('ERROR: ') + 'The name "noise" is not available in extras.\n'
            msg += red('     : ') + 'We need the RV photon noise to scale the noise in the diagnoses.'
            clogger.fatal(msg)
            return       

    def transform_to_sane_units(self):
        """
        Transform the measurements to units that can be plotted easily,
        i.e. are of order unity (or close). This changes RVs, BIS.
        FWHM in km/s is closer to unity than in m/s.
        Returns (weighted_average_RV, weighted_average_BIS)
        """
        self.set_errors_indicators()

        if self.units == 'km/s':
            weighted_average_RV = np.average(self.vrad, weights=1./self.error**2)
            weighted_average_BIS = np.average(self.extras.bis_span, weights=1./self.error**2)

            # subtract average and transform to m/s
            self.vrad = (self.vrad - weighted_average_RV) * 1e3
            self.error *= 1e3
            self.units = 'm/s'

            # convert the BIS to m/s
            newBIS = (self.extras.bis_span - weighted_average_BIS) * 1e3
            newBISsig = self.extras.sig_bis * 1e3
            self.extras = self.extras._replace(bis_span = newBIS, sig_bis=newBISsig)


            msg = blue('INFO: ') + 'Subtracted weighted average (%f km/s) from RV\n' % weighted_average_RV
            msg += blue('    : ') + '                            (%f km/s) from BIS\n' % weighted_average_BIS
            msg += blue('INFO: ') + 'Converted to m/s.'
            clogger.info(msg)


        return weighted_average_RV, weighted_average_BIS


    def save_before_after_fibers(self):
        """ Save the observations in two separate files before and after fiber change. """
        
        time_of_offset = 57170  # HARPS upgrade
        assert (self.time > time_of_offset).any(), 'No observations after BJD=%d' % time_of_offset

        # the index of the first observation after time_of_offset
        ind_change = np.where(self.time > time_of_offset)[0][0]

        # assume only one file
        assert len(self.provenance) == 1
        original_filename = os.path.basename(self.provenance.keys()[0])
        filename1 = original_filename.split('_')
        filename1.insert(1, '_BeforeOffset_')
        filename1 = ''.join(filename1)

        msg = blue('INFO: ') + 'Writing observations before %d to file "%s"' % (time_of_offset, filename1)
        clogger.info(msg)
        #
        rvIO.write_rv(self, filename1, imax=ind_change)

        filename2 = original_filename.split('_')
        filename2.insert(1, '_AfterOffset_')
        filename2 = ''.join(filename2)

        msg = blue('INFO: ') + 'Writing observations after %d to file "%s"' % (time_of_offset, filename2)
        clogger.info(msg)
        #
        rvIO.write_rv(self, filename2, imin=ind_change)


    def do_plot_obs(self, newFig=True, axis=None, leg=True, save=None, offsets=None, LP=False):
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
        ax2.plot(years, self.vrad.mean() * np.ones_like(years), alpha=0) # Create a dummy plot

        # plot each files' values
        if offsets:
            assert isinstance(offsets, list)
            assert len(offsets) == len(self.provenance)

        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
            #                     height=err[:m], beam=0.5,\
            #                     pen=pg.mkPen(None))
                                # pen={'color': 0.8, 'width': 2})
            # p.addItem(e)
            # p.plot(t[:m], rv[:m], symbol='o')
            offs = offsets[i] if offsets else 0.
            ax1.errorbar(t[:m], rv[:m]+offs, yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err = t[m:], rv[m:], err[m:]
        
        if LP:
            # this is 1 October 2012
            ax1.axvline(x=56202, ls='--', color='k')

        if leg: ax1.legend()
        plt.tight_layout()
        ax2.ticklabel_format(useOffset=False)

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)

        return fig
        # plt.show()
        # pg.QtGui.QApplication.exec_()

    def do_plot_obs_pretty(self, newFig=True, leg=False, save=None, offsets=None, show_years=False):
        """ Plot the observed radial velocities as a function of time.
        Data from each file are color coded and labeled. Pretty plot, ready for publication.
        """
        full_path = self.provenance.keys()[0]
        bn = os.path.basename(full_path)
        i = bn.rfind('_harps_mean_corr.rdb')
        if i == -1:
            i = bn.rfind('_harps_mean.rdb')
        if i == -1:
            i = bn.rfind('_harps.rdb')
        star = bn[:i]

        with plt.rc_context(params):

            class MyFormatter(ticker.ScalarFormatter):
                def __call__(self, x, pos=None):
                    # call the original LogFormatter
                    rv = ticker.ScalarFormatter.__call__(self, x, pos)
                    # check if we really use TeX
                    if plt.rcParams["text.usetex"]:
                        # if we have the string ^{- there is a negative exponent
                        # where the minus sign is replaced by the short hyphen
                        rv = re.sub(r'-', r'\mhyphen', rv)
                    return rv


            figwidth = 3.543311946  # in inches = \hsize = 256.0748pt
            figheight = 0.75 * figwidth
            # this is Seaborn's "colorblind" pallete
            # colors = ['#0072b2', '#009e73', '#d55e00', '#cc79a7', '#f0e442', '#56b4e9']
            # this is Seaborn's "muted" pallete
            colors = ['#4878cf', '#6acc65', '#d65f5f', '#b47cc7', '#c4ad66', '#77bedb']
            t, rv, err = self.time-50000, self.vrad, self.error # temporaries
            
            if newFig: 
                fig = plt.figure(figsize=(figwidth, figheight))

            ax2 = fig.add_subplot(111)
            ax2.set_title(star, loc='right', fontsize=params['font.size'])
            lpad = 20 if show_years else 8
            ax2.set_xlabel('BJD - 2450000 [days]', labelpad=lpad)
            ax2.set_ylabel(r'RV [$\ms$]')
            ax1 = ax2
            
            if show_years:
                ny, years = self.get_years_observations()
                ax1 = ax2.twiny()
                ax1.xaxis.tick_bottom()
                ax2.xaxis.tick_top()
                ax2.plot(years, np.ones_like(years), alpha=0) # Create a dummy plot

            # plot each files' values
            if offsets:
                assert isinstance(offsets, list)
                assert len(offsets) == len(self.provenance)

            for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
                m = n-nout # how many values are there after restriction
                offs = offsets[i] if offsets else 0.

                ax1.errorbar(t[:m], rv[:m]+offs, yerr=err[:m],
                             fmt='o', color=colors[i], 
                             mec='none', ms=2, capsize=0, elinewidth=0.5)
                t, rv, err = t[m:], rv[m:], err[m:]
            
            plt.tight_layout()
            ax2.ticklabel_format(useOffset=False)
            ax1.yaxis.set_major_formatter(MyFormatter())
            print ax1.margins()
            ax1.margins(0.05)

            if save:
                msg = yellow('INFO: ') + 'Saving figure to %s' % save
                clogger.info(msg)
                plt.savefig(save, bbox_inches='tight')

        ax2.set_ylabel('RV [%s]' % self.units) # otherwise, matplotlib complains...
        plt.show()

        return fig
        # plt.show()
        # pg.QtGui.QApplication.exec_()


    def do_plot_obs_together(self, q=None, newFig=True, leg=False, save=None, offsets=None, LP=False):
        """ Plot the observed radial velocities as a function of time, together with other variables.
        Data from each file are color coded and labeled.
        """

        ## handle inexistent field
        q = q.split(',')
        print q
        nvar = len(q)
        for qq in q:
            print qq
            if qq not in self.extras._fields:
              from shell_colors import red
              msg = red('ERROR: ') + 'The name "%s" is not available in extras.\n' % qq
              clogger.fatal(msg)
              return

        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = self.time, self.vrad, self.error # temporaries
        
        if newFig: 
            fig = plt.figure(figsize=(8,5+3*nvar))
        ax2 = fig.add_subplot(1+nvar,1,1)
        ax2.set_ylabel('RV [%s]'%self.units)

        ny, years = self.get_years_observations()
        ax1 = ax2.twiny()
        ax1.xaxis.tick_bottom()
        ax2.xaxis.tick_top()
        ax2.plot(years, self.vrad.mean() * np.ones_like(years), alpha=0) # Create a dummy plot

        # plot each files' values
        if offsets:
            assert isinstance(offsets, list)
            assert len(offsets) == len(self.provenance)

        for i, (fname, [n, nout]) in enumerate(sorted(self.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # e = pg.ErrorBarItem(x=t[:m], y=rv[:m], \
            #                     height=err[:m], beam=0.5,\
            #                     pen=pg.mkPen(None))
                                # pen={'color': 0.8, 'width': 2})
            # p.addItem(e)
            # p.plot(t[:m], rv[:m], symbol='o')
            offs = offsets[i] if offsets else 0.
            ax1.errorbar(t[:m], rv[:m]+offs, yerr=err[:m], \
                         fmt='o'+colors[i], label=fname)
            t, rv, err = t[m:], rv[m:], err[m:]
        
        if LP:
            # this is 1 October 2012
            ax1.axvline(x=56202, ls='--', color='k')

        if leg: ax1.legend()
        ax2.ticklabel_format(useOffset=False)

        ## plot the other quantities
        for i, qq in enumerate(q):
            ax3 = fig.add_subplot(1+nvar, 1, 2+i, sharex=ax1)
            if qq=='rhk': 
                ax3.errorbar(self.time, self.extras.rhk, yerr=self.extras.sig_rhk, fmt='o')
            elif qq=='fwhm':
                f = 2.35e-3 if self.units=='m/s' else 2.35
                ax3.errorbar(self.time, self.extras.fwhm, yerr=f*self.error, fmt='o')
            else:
                ind = self.extras._fields.index(qq) # index corresponding to this quantity
                plt.plot(self.time, self.extras[ind], 'o', label=qq)

            ax3.set_ylabel(qq)

        ax3.set_xlabel('Time [days]', labelpad=20)
        plt.tight_layout()
        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)

        return fig
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
        # elif extra == 'fwhm' and ask_yes_no('Should I remove a linear trend first? (y/N) ', False):
        #     m, b = np.polyfit(t, self.extras[i], 1)
        #     yp = np.polyval([m, b], t)
        #     plt.plot(t, self.extras[i]-yp, 'o', label=extra+' (linear trend removed)')
        elif extra == 'fwhm':
            f = 2.35e-3 if self.units=='m/s' else 2.35
            plt.errorbar(t, self.extras.fwhm, yerr= f*self.error, fmt='o', label=extra)
        else:
            plt.plot(t, self.extras[i], 'o', label=extra)
        
        plt.xlabel('Time [days]')
        plt.ylabel(extra + ' []')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05))
        plt.minorticks_on()
        plt.tight_layout()

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

        par = self.fit['params']
        keplerians = int(len(par)/6)
        P, K, ecc, omega, T0, gam = [par[j::6] for j in range(6)]
        gam = gam[0]

        tt = np.linspace(self.time.min(), self.time.max(), 100*np.ptp(self.time)/P)
        final = np.zeros_like(tt)

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

        # plt.show()

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
            t, r, err = t[m:], r[m:], err[m:]
        
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

    def get_time_to_plot(self, P=None, oversample=10):
        """
        Returns sufficiently resolved time vector for plots
        """
        std = self.time.std()
        N = len(self.time)
        minim = self.time.min() - 1.*std
        maxim = self.time.max() + 1.*std
        if P is None:
            return np.linspace(minim, maxim, oversample*N)
        else:
            return np.linspace(minim, maxim, oversample*(maxim-minim)/P)

    def is_in_extras(self, extra):
        return extra in self.extras._fields


class BasicTimeSeries():
    time = None
    vrad = None
    error = None


perc01 = 0.001 # 0.1% FAP
perc1 = 1.  # 1% FAP
perc10 = 10.  # 10% FAP

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
    star_name = None # name of the star


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
        from tqdm import tqdm
        from time import time
        from multiprocessing import cpu_count

        ncpu = cpu_count()
        name = '_' + self.__class__.__name__

        if name == '_gls':
            from OPEN.ext.glombscargle import glombscargle
        else:
            return

        # temporaries
        temp_per = copy(self)

        f = temp_per.freq
        omegas = 2.*np.pi*temp_per.freq
        p = temp_per._upow

        perm = 1000 # int(1000/perc1) # (use 1000 for 1% fap or 10000 for 0.1% fap)

        try:
            peaks = self.peaks
            if len(self.peaks) != perm: raise AttributeError
        except AttributeError:
            t1 = time()
            maxPowers = []
            t = temp_per.t
            for k in tqdm(xrange(perm)):
                permutted = np.random.permutation(zip(temp_per.y, temp_per.error))
                # temp_per.y = permutted[:,0]
                y = permutted[:,0]
                # temp_per.error = permutted[:,1]
                err = permutted[:,1]

                power = glombscargle(t, y, err, omegas, ncpu)[0]
                # calc()
                # temp_per._plot_pg()
                # powermaxP = temp_per.power.max()
                powermaxP = power.max()
                maxPowers.append(powermaxP)
        
            peaks = np.sort(maxPowers)

        return peaks


    def _plot(self, doFAP=False, dobFAP=False, faps=None, verts=None, newFig=True, axes=None, save=None, **kwargs):
        """
        Plot this periodogram.
        """
        xlabel = 'Period [d]'
        ylabel = 'Power'
        do_title = kwargs.pop('title', True)
        do_legend = kwargs.pop('legend', True)
        do_labels = kwargs.pop('labels', True)
        color = kwargs.pop('color', 'b')
        FAPcolor = kwargs.pop('FAPcolor', 'k')

        if newFig and save:
            self.fig = plt.figure(figsize=(8,4))
        elif newFig:
            self.fig = plt.figure()
        else:
            self.fig = plt.gcf()

        if axes is None:
            self.ax = self.fig.add_subplot(1,1,1)
        else:
            self.ax = axes

        if do_title:
            self.ax.set_title("Normalized periodogram")

        if do_labels:
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel(ylabel)
        if self.power.max() < 1e-300:  # apparently, Metplotlib can't plot these small values
            clogger.warning(yellow('Warning: ')+'Max value < 1e-300, plotting normalized periodogram')
            self.ax.semilogx(1./self.freq, self.power/self.power.max(), color=color, ls='-', **kwargs)
        else:
            self.ax.semilogx(1./self.freq, self.power, color=color, ls='-', **kwargs)
        # plot FAPs
        if doFAP:
            # do default FAPs of 10%, 1% and 0.1%
            if faps is None: clogger.warning(yellow('Warning: ')+'Plotting default FAPs')

            pmin = 1./self.freq.min()
            pmax = 1./self.freq.max()
            plvl1 = self.powerLevel(0.1) # 10% FAP
            plvl2 = self.powerLevel(0.01) # 1% FAP
            plvl3 = self.powerLevel(0.001) # 0.1% FAP
            self.ax.axhline(y=plvl1, color=FAPcolor, ls='-', label='10%')
            self.ax.axhline(y=plvl2, color=FAPcolor, ls='--', label='1%')
            self.ax.axhline(y=plvl3, color=FAPcolor, ls=':', label='0.1%')
            if do_legend: 
                self.ax.legend(frameon=True)
        if dobFAP:
            # calculate FAP by bootstrap
            self.peaks = self.FAP_by_bootstrap()

            index01 = int( ((1-perc01/100.0) * len(self.peaks)) )
            index1 = int( ((1-perc1/100.0) * len(self.peaks)) )
            index10 = int( ((1-perc10/100.0) * len(self.peaks)) )
            self.powerFAP_01 = self._normalize_value(self.peaks[index01])
            self.powerFAP_1 = self._normalize_value(self.peaks[index1])
            self.powerFAP_10 = self._normalize_value(self.peaks[index10])

            clogger.info(blue('INFO: ')+'Plotting bootstrap FAPs')
    
            self.ax.axhline(self.powerFAP_10, c=FAPcolor, lw=2, ls='-', label='10%')
            self.ax.axhline(self.powerFAP_1, c=FAPcolor, lw=2, ls='--', label='1%')
            self.ax.axhline(self.powerFAP_01, c=FAPcolor, lw=2, ls=':', label='0.1%')

            if do_legend: 
                self.ax.legend(frameon=True)

        # plot vertical lines
        if verts is not None:
            for v in verts:
                self.ax.axvline(x=v, color='k', ls='--', lw=2, alpha=0.5, label="_nolegend_") 

        # plt.tight_layout()

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)


    def _plot_pretty(self, doFAP=False, dobFAP=False, faps=None, verts=None, save=None, **kwargs):
        """ Plot this periodogram. Ready for publication. """

        # this is Seaborn's "colorblind" pallete
        colors = ['#0072b2', '#009e73', '#d55e00', '#cc79a7', '#f0e442', '#56b4e9']
        # this is Seaborn's "muted" pallete
        # colors = ['#4878cf', '#6acc65', '#d65f5f', '#b47cc7', '#c4ad66', '#77bedb']

        xlabel = 'Period [days]'
        ylabel = 'Power'
        do_title = kwargs.pop('title', True)
        do_labels = kwargs.pop('labels', True)
        color = kwargs.pop('color', colors[0])
        FAPcolor = kwargs.pop('FAPcolor', 'k')

        with plt.rc_context(params):
            figwidth = 3.543311946  # in inches = \hsize = 256.0748pt
            figheight = 0.5 * figwidth


            fig = plt.figure(figsize=(figwidth, figheight))
            ax = fig.add_subplot(111)
            ax.set_title(self.star_name, loc='right', fontsize=params['font.size'])

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if self.power.max() < 1e-300:  # apparently, Metplotlib can't plot these small values
                clogger.warning(yellow('Warning:')+' Max value < 1e-300, plotting normalized periodogram')
                ax.semilogx(1./self.freq, self.power/self.power.max(), color=color, ls='-', **kwargs)
            else:
                ax.semilogx(1./self.freq, self.power, color=color, ls='-', lw=0.5, **kwargs)
            # plot FAPs
            if doFAP:
                # do default FAPs of 10%, 1% and 0.1%
                if faps is None: clogger.warning(yellow('Warning: ')+'Plotting default FAPs')

                pmin = 1./self.freq.min()
                pmax = 1./self.freq.max()
                plvl1 = self.powerLevel(0.1) # 10% FAP
                plvl2 = self.powerLevel(0.01) # 1% FAP
                plvl3 = self.powerLevel(0.001) # 0.1% FAP
                ax.axhline(y=plvl1, color='k', ls='-', label='10%')
                ax.axhline(y=plvl2, color='k', ls='--', label='1%')
                ax.axhline(y=plvl3, color='k', ls=':', label='0.1%')
                if do_legend: 
                    ax.legend(frameon=True)
            if dobFAP:
                # calculate FAP by bootstrap
                self.peaks = self.FAP_by_bootstrap()

                index01 = int( ((1-perc01/100.0) * len(self.peaks)) )
                index1 = int( ((1-perc1/100.0) * len(self.peaks)) )
                index10 = int( ((1-perc10/100.0) * len(self.peaks)) )
                self.powerFAP_01 = self._normalize_value(self.peaks[index01])
                self.powerFAP_1 = self._normalize_value(self.peaks[index1])
                self.powerFAP_10 = self._normalize_value(self.peaks[index10])

                clogger.info(blue('INFO: ')+'Plotting bootstrap FAPs (1% and 0.1%)')
        
                # ax.axhline(self.powerFAP_10, xmin=0.05, xmax=0.95, c=FAPcolor, lw=0.5, ls='-', label='10%')
                ax.axhline(self.powerFAP_1, xmin=0.04, xmax=0.96, c=FAPcolor, lw=0.5, ls='--', dashes=(4,4),  label='1%')
                ax.axhline(self.powerFAP_01, xmin=0.04, xmax=0.96, c=FAPcolor, lw=0.5, ls=':', label='0.1%')


            # plot vertical lines
            if verts is not None:
                for v in verts:
                    ax.axvline(x=v, color='k', ls='--', lw=2, alpha=0.5, label="_nolegend_") 

            try:
                if ax.get_ylim()[1] > self.powerFAP_01+2:
                    pass
                else:
                    ax.set_ylim([0, self.powerFAP_01+2])
            except:
                pass
            fig.tight_layout()

            if save:
                msg = yellow('INFO: ') + 'Saving figure to %s' % save
                clogger.info(msg)
                plt.savefig(save)
            plt.show()



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
      bbin = np.argmax(self.power)
      # Maximum power
      pmax = self._upow[bbin]

      # rms = np.sqrt(self._YY * (1.-pmax))

      # Get the curvature in the power peak by fitting a parabola y=aa*x^2
      if (bbin > 1) and (bbin < len(self.freq)-2):
        # Shift the parabola origin to power peak
        xh = (self.freq[bbin-1:bbin+2] - self.freq[bbin])**2
        yh = self._upow[bbin-1:bbin+2] - self._upow[bbin]
        # Calculate the curvature (final equation from least square)
        aa = np.sum(yh*xh)/np.sum(xh*xh)
        nt = float(self.N)
        f_err = np.sqrt(-2./nt * pmax/aa*(1.-pmax)/pmax)
        Psin_err = np.sqrt(-2./nt* pmax/aa*(1.-pmax)/pmax) / self.freq[bbin]**2
      else:
        f_err = None
        Psin_err = None

      fbest = self.freq[bbin]
      # amp = np.sqrt(self._a[bbin]**2 + self._b[bbin]**2)
      # ph  = np.arctan2(self._a[bbin], self._b[bbin]) / (2.*pi)
      # T0  = np.min(self.th) - ph/fbest
      # Re-add the mean
      # offset = self._off[bbin] + self._Y

      # Statistics
      print "Generalized LS - statistical output"
      print 33*"-"
      if verbose:
        print "Number of input points:     %6d" % (nt)
        # print "Weighted mean of dataset:   % e" % (self._Y)
        # print "Weighted rms of dataset:    % e" % (sqrt(self._YY))
        print "Time base:                  % e" % (max(self.th) - min(self.th))
        print "Number of frequency points: %6d" % (len(self.freq))
        print
        print "Maximum power, p :    % e " % (self.power[bbin])
        print "Maximum power (without normalization):   %e" % (pmax)
        print "Normalization    : ", self.norm
        # print "RMS of residuals :    % e " % (rms)
        if self.error is not None:
          print "  Mean weighted internal error:  % e" %(np.sqrt(nt/np.sum(1./self.error**2)))
        print "Best sine frequency : % e +/- % e" % (fbest, f_err)
        print "Best sine period    : % e +/- % e" % (1./fbest, Psin_err)
        # print "Amplitude:          : % e +/- % e" % (amp, sqrt(2./nt)*rms)
        # print "Phase (ph) : % e +/- % e" % (ph, sqrt(2./nt)*rms/amp/(2.*pi))
        # print "Phase (T0) : % e +/- % e" % (T0, sqrt(2./nt)*rms/amp/(2.*pi)/fbest)
        # print "Offset     : % e +/- % e" % (offset, sqrt(1./nt)*rms)
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
        for i, planet in enumerate(list(ascii_lowercase)[1:self.nplanets+1]):
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
        # plt.show()

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
        # plt.show()

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

        if self.context[2] == '2':  # model with jitter
            self.jitter = True
        else:
            self.jitter = False

        # number of planets (this is set in the namelist so it should be bullet-proof)
        self.nplanets = int(context[1])
        # number of obervatories (we avoid reading the input.rv file)
        self.nobserv = self.npar - 5*self.nplanets # subtract planets' parameters
        if self.jitter: self.nobserv -= 1 # subtract jitter parameter
        if self.gp: self.nobserv -= 5 # subtract hyperparameters (assumed 5 for now)

        self.read_stats_file()

    def read_stats_file(self):
        """ 
        Read stats.dat file for mean, ML and MAP parameters
        """
        with open(self.root+'stats.dat') as f:
            stats = f.readlines()

        npar = self.npar
        # if self.gp:
        #     self.nplanets = (npar-4)/5
        # elif self.gp and (npar == 5): 
        #     self.nplanets = 0
        #     self.gp_only = True
        # else:
        #     self.nplanets = npar/5
        

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

        if self.jitter:
            i = -self.nobserv-1
            self.jitter_map = self.par_map.pop(i)
            self.jitter_mle = self.par_mle.pop(i)
            self.jitter_mean = self.par_mean.pop(i)
            self.jitter_sigma = self.par_sigma.pop(i)


        self.only_vsys = False
        j = 1 if self.jitter else 0
        j = j+5 if self.gp else j
        if self.npar-j == self.nobserv:  # systematic velocity only
            self.only_vsys = True

    def read_gp_file(self):
        filename = self.root+'gp.dat'
        try:
            self.pred_t, self.pred_y, self.pred_std = np.loadtxt(filename, unpack=True)
        except IOError:
            self.pred_t, self.pred_y, self.pred_std = (np.nan, np.nan, np.nan)

    def read_iterations(self):
        """
        Read the [root]ev.dat and [root]phys_live.points files which contain
        the set of rejected points and the current set of live points, respectively.
        First npar columns are the parameter values for each iteration + live ones
        npar+1 column are the log-likelihood values
        """
        filename1 = self.root + 'ev.dat'
        filename2 = self.root + 'phys_live.points'

        rejected_points = np.genfromtxt(filename1, unpack=True)
        live_points = np.genfromtxt(filename2, unpack=True)

        j = 1 if self.jitter else 0
        j = j+5 if self.gp else j

        self.posterior_samples = np.append(rejected_points[:self.npar, :], live_points[:self.npar, :], axis=1)

        parameter_names = ['period', 'K', 'ecc', 'omega', 't0']
        parameter_names = np.array([[par + '_'+str(i) for par in parameter_names] for i in range(self.nplanets)])

        # we have to name every namedtuple class differently
        p = namedtuple('parameter'+str(self.nplanets), parameter_names.flat)
        # and assign it to the module globals, so that pickle doesn't complain
        globals()[p.__name__] = p

        # this makes a namedtuple of all the parameter values
        self.trace = p._make(self.posterior_samples[:-1-j,:])


    def get_fit_statistics(self, system):
        """
        Read the [root]phys_live.points file which contains the current set of
        live points, to obtain the maximum likelihood value.
        """
        filename = self.root + 'phys_live.points'
        live_points = np.genfromtxt(filename, unpack=True)

        max_llike = max(live_points[-2])
        k = live_points.shape[0]
        n = system.time.size

        self.AIC = 2*k - 2*max_llike
        self.BIC = k*np.log(n) -2*max_llike
        self.logL = max_llike

    def compress_chains(self):
        """
        Save all MultiNest output files in a zip file, allowing for future restarts.
        Returns the zip file name.
        """
        tstr = strftime("%Y%m%d-%H%M%S")
        zfilename = 'chains/nest-out-'+tstr+'.zip'
        nest_output_files = set(glob.glob(self.root+'*')) - set(glob.glob(self.root+'*.zip'))

        with ZipFile(zfilename, 'a', compression=8) as zf:
            for f in nest_output_files:
                zf.write(f)
            # include the namelist used to run MultiNest
            namelist_file = os.path.join(os.path.dirname(__file__), 'multinest/namelist1')
            zf.write(namelist_file)

        msg = blue('INFO: ') + 'Saved output files to %s' % zfilename
        clogger.info(msg)

        return zfilename

    def print_best_solution(self, system):

        par_mean, par_sigma, par_mle, par_map = self.par_mean, self.par_sigma, self.par_mle, self.par_map

        try:
            self.AIC
        except AttributeError:
            self.get_fit_statistics(system)


        msg = yellow('RESULT: ') + 'Parameters summary'
        clogger.info(msg)

        print '%8s %14s %9s %14s %14s' % ('', 'mean', '+- sigma', 'ML', 'MAP')
        ## loop over planets
        i = -1
        for i, planet in enumerate(list(ascii_lowercase)[1:self.nplanets+1]):
            print yellow(planet)
            print '%8s %14.3f %9.3f %14.3f %14.3f' % ('P',     par_mean[5*i], par_sigma[5*i], par_mle[5*i], par_map[5*i])
            print '%8s %14.3f %9.3f %14.3f %14.3f' % ('K',     par_mean[5*i+1], par_sigma[5*i+1], par_mle[5*i+1], par_map[5*i+1])
            print '%8s %14.3f %9.3f %14.3f %14.3f' % ('ecc',   par_mean[5*i+2], par_sigma[5*i+2], par_mle[5*i+2], par_map[5*i+2])
            print '%8s %14.3f %9.3f %14.3f %14.3f' % ('omega', par_mean[5*i+3], par_sigma[5*i+3], par_mle[5*i+3], par_map[5*i+3])
            print '%8s %14.2f %9.2f %14.2f %14.2f' % ('t0',    par_mean[5*i+4], par_sigma[5*i+4], par_mle[5*i+4], par_map[5*i+4])
        
            print 

            from .utils import mjup2mearth, mearth2msun, mean_sidereal_day, au2m
            P, K, ecc = par_map[5*i], par_map[5*i+1], par_map[5*i+2]
            P_error, K_error, ecc_error = par_sigma[5*i], par_sigma[5*i+1], par_sigma[5*i+2]

            try:
                from uncertainties import ufloat
                from uncertainties.umath import sqrt
                P = ufloat(P, P_error)
                K = ufloat(K, K_error)
                ecc = ufloat(ecc, ecc_error)
                m_mj = 4.919e-3 * system.star_mass**(2./3) * P**(1./3) * K * sqrt(1-ecc**2)
                m_me = m_mj * mjup2mearth
                a = ((system.star_mass + m_me*mearth2msun)/(m_me*mearth2msun)) * sqrt(1.-ecc**2) * K * (P*mean_sidereal_day/(2*np.pi)) / au2m

                print '%8s %11.4f +- %7.4f [MJup]  %11.4f +- %7.4f [MEarth]' % ('m sini', m_mj.n, m_mj.s, m_me.n, m_me.s)
                print '%8s %11.4f +- %7.4f [AU]' % ('a', a.n, a.s)

            except ImportError:
                m_mj = 4.919e-3 * system.star_mass**(2./3) * P**(1./3) * K * np.sqrt(1-ecc**2)
                m_me = m_mj * mjup2mearth
                a = ((system.star_mass + m_me*mearth2msun)/(m_me*mearth2msun)) * np.sqrt(1.-ecc**2) * K * (P*mean_sidereal_day/(2*np.pi)) / au2m

                print '%8s %11.3f [MJup] %11.3f [MEarth]' % ('m sini', m_mj, m_me)
                print '%8s %11.3f [AU]' % ('a', a)

            print 

        print yellow('system')
        if self.jitter:
            # jitter parameter
            print '%8s %14.3f %9.3f %14.3f %14.3f' % ('jitter', self.jitter_mean, self.jitter_sigma, self.jitter_mle, self.jitter_map)
        if self.gp:
            # in this case, the vsys is before the hyperparameters
            print '%8s %14.3f %9.3f %14.3f %14.3f' % ('vsys', par_mean[5*i+5], par_sigma[5*i+5], par_mle[5*i+5], par_map[5*i+5])
            print yellow('GP')
            for j in sorted(range(1, 6), reverse=True):
                print '%8s %14.3f %9.3f %14.3f %14.3f' % ('sigma'+str(6-j), par_mean[-j], par_sigma[-j], par_mle[-j], par_map[-j])
        else:
            # in this case, the vsys parameters are the last ones
            nobs = self.nobserv
            for i in range(nobs):
                j = -nobs+i
                print '%8s %14.3f %9.3f %14.3f %14.3f' % ('offs'+str(i+1), par_mean[j], par_sigma[j], par_mle[j], par_map[j])

        print yellow('\nfit')
        print '  logL: %8.3f   AIC: %8.3f   BIC:%8.3f' % (self.logL, self.AIC, self.BIC)
        print

    def confidence_intervals(self):
        try:
            self.trace
        except AttributeError:
            self.read_iterations()

        for i, name in enumerate(self.trace._fields):
            print '%10s' % name,
            print '  ', 'mean:', np.mean(self.trace[i]), 'std:', np.std(self.trace[i])
            if (i+1)%5==0: print

    def save_fit_to(self, system):
        
        if system.fit is not None: # there is already a fit
            clogger.warning(yellow('Warning: ')+'Replacing older fit')

        system.fit = {}
        
        t, rv, err = system.time, system.vrad, system.error # temporaries
        vel = np.zeros_like(t)
        
        # if only one observatory, observ is always 1
        observ = np.ones_like(vel, dtype=int)
        # else it takes a different integer value for each observatory
        if self.nobserv > 1:
            n_each_observ = [v[0] for v in system.provenance.values()]
            chunks = np.split(observ, np.cumsum(n_each_observ))[:-1]
            for i in range(self.nobserv):
                chunks[i] = chunks[i] * (i+1)

            observ = np.concatenate(chunks)

        j = 5 if self.gp else 0
        if self.nobserv > 1:
            vsys = self.par_map[-self.nobserv:]
        else:
            vsys = [self.par_map[-1-j]]


        ## MAP estimate of the parameters
        if self.gp and not self.only_vsys:
            par_map = self.par_map[:-5] 
            hyper_map = self.par_map[-5:]
            if self.gp_only:
                vel = gp_predictor(t, rv, err, par_map, hyper_map, 'constant')
            else:
                vel = gp_predictor(t, rv, err, par_map, hyper_map, 'keplerian')
        elif not self.only_vsys:
            par_map = self.par_map

            # only planet(s)' parameters
            planets_par_map = par_map[:-self.nobserv]

            P = planets_par_map[::5]
            K = planets_par_map[1::5]
            ecc = planets_par_map[2::5]
            omega = planets_par_map[3::5]
            t0 = planets_par_map[4::5]
            # vsys goes as 0. into get_rvn
            par = [P, K, ecc, omega, t0, 0.]

            args = [t] + par + [vel]
            get_rvn(*args)

        # we add the velocity offsets here
        for i in range(self.nobserv):
            ind = np.where(observ==(i+1))[0]
            vel[ind] += vsys[i]

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
        # plt.show()

    def do_plot_trace(self, parameter):
        try:
            self.trace
        except AttributeError:
            self.read_iterations()

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

    def do_plot_map(self, system, legend=True, save=None, oversample=10):
        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = system.time, system.vrad, system.error # temporaries
        # tt = system.get_time_to_plot(P=self.par_map[0])
        tt = system.get_time_to_plot(oversample=oversample)
        vel = np.zeros_like(tt)
        velt = np.zeros_like(t)

        # if only one observatory, observ is always 1
        observ = np.ones_like(velt, dtype=int)
        # else it takes a different integer value for each observatory
        if self.nobserv > 1:
            n_each_observ = [v[0] for v in system.provenance.values()]
            chunks = np.split(observ, np.cumsum(n_each_observ))[:-1]
            for i in range(self.nobserv):
                chunks[i] = chunks[i] * (i+1)

            observ = np.concatenate(chunks)

        # get the RV offsets
        j = 5 if self.gp else 0
        if self.nobserv > 1:
            vsys = self.par_map[-self.nobserv:]
        else:
            vsys = [self.par_map[-1-j]]

        # we add the velocity offsets here
        # residuals will be RV-velt
        for i in range(self.nobserv):
            # ind = np.where(observ==(i+1))[0]
            # vel[ind] += vsys[i]
            ind = np.where(observ==(i+1))[0]
            velt[ind] += vsys[i]

        ## MAP estimate of the parameters
        if self.gp:
            par_map = self.par_map[:-5]
            hyper_map = self.par_map[-5:]
            # print par_map, hyper_map
            pred = gp_predictor(t, rv, err, par_map, hyper_map, 'stub')
        else:
            par_map = self.par_map

        newFig=True
        if newFig:
            fig = plt.figure()

        gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharex=ax1)

        if self.gp and not self.only_vsys:
            # plot best solution
            # only planet(s)' parameters
            planets_par_map = self.par_map[:-5-self.nobserv]

            P = planets_par_map[::5]
            K = planets_par_map[1::5]
            ecc = planets_par_map[2::5]
            omega = planets_par_map[3::5]
            t0 = planets_par_map[4::5]
            par = [P, K, ecc, omega, t0, 0.]

            args = [tt] + par + [vel]
            get_rvn(*args)
            ax1.plot(tt, vel+vsys, '-g', lw=2.5, label='MAP')            
        elif self.only_vsys:
            pass
        else:
            # plot best solution
            # only planet(s)' parameters
            planets_par_map = par_map[:-self.nobserv]

            P = planets_par_map[::5]
            K = planets_par_map[1::5]
            ecc = planets_par_map[2::5]
            omega = planets_par_map[3::5]
            t0 = planets_par_map[4::5]
            par = [P, K, ecc, omega, t0, 0.]

            args = [tt] + par + [vel]
            get_rvn(*args)
            ax1.plot(tt, vel+vsys[0], '-g', lw=2.5, label='MAP')
            # for t in t0:
            # t0s = t0 + P * np.arange(-150, 51)
            # tts = get_tt(P[0], ecc[0], omega[0], t0s)
            # ax1.vlines(x=t0s, ymin=-500, ymax=500, color='k', alpha=0.8)
            # ax1.vlines(x=tts, ymin=-500, ymax=500, color='r', alpha=0.8)

        # plot GP predictions
        if self.gp:
            # ax1.plot(t, pred, 'ro')
            pred_y = self.pred_y+vsys[0]
            ax1.plot(self.pred_t, pred_y, '-k', lw=1.5, label='GP mean')
            ax1.fill_between(self.pred_t, y1=pred_y-2*self.pred_std, y2=pred_y+2*self.pred_std,
                                          color='k', alpha=0.3, label='2*std')

        # vel = np.zeros_like(t)
        if self.gp or self.only_vsys:
            # vel = self.par_map[0] * np.ones_like(tt)
            pass
        else:
            args = [t] + par + [velt]
            get_rvn(*args)

        for i, (fname, [n, nout]) in enumerate(sorted(system.provenance.iteritems())):
            m = n-nout # how many values are there after restriction
            
            # plot each files' values offset by systematic velocities
            if self.gp:
                # here we don't remove vsys yet but this needs to be taken care of!!
                ax1.errorbar(t[:m], rv[:m], yerr=err[:m], fmt='o'+colors[i], label=os.path.basename(fname))
            else:
                # because we add each offset to the model RV, we subtract them from the observations
                ax1.errorbar(t[:m], rv[:m], yerr=err[:m], fmt='o'+colors[i], label=os.path.basename(fname))
            # plot residuals
            if self.gp:
                ax2.errorbar(t[:m], rv[:m]-pred[:m], yerr=err[:m], fmt='o'+colors[i], label=fname)
                pred = pred[m:]
            else:
                ax2.errorbar(t[:m], rv[:m]-velt[:m], yerr=err[:m], fmt='o'+colors[i], label=fname)
                velt = velt[m:]

            t, rv, err = t[m:], rv[m:], err[m:]



        # plot systematic velocity
        for i, v in enumerate(vsys):
            ax1.axhline(y=v, ls='--', color=colors[i], alpha=0.3)
        ax2.axhline(y=0., ls='--', color='k', alpha=0.3)

        ax2.set_xlabel('Time [days]')
        ax1.set_ylabel('RV [%s]'%system.units)
        ax2.set_ylabel('Residuals [%s]'%system.units)
        if legend: ax1.legend()
        fig.tight_layout()
        # fig.ticklabel_format(useOffset=False)

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            fig.savefig(save)


    def do_plot_map_phased(self, system, legend=False, plot_gp=True, save=None, oversample=10):
        # if systematic velocity only, there is nothing to do here
        if self.only_vsys: return

        colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = system.time, system.vrad, system.error # temporaries
        tt = system.get_time_to_plot(oversample=oversample)
        vel = np.zeros_like(tt)

        # get the RV offsets
        j = 4 if self.gp else 0
        if self.nobserv > 1:
            vsys = self.par_map[-self.nobserv:]
        else:
            vsys = [self.par_map[-1-j]]

        ## MAP estimate of the parameters
        if self.gp:
            par_map = self.par_map[:-4]  # don't care about the hyperparameters for now
        else:
            par_map = self.par_map

        newFig=True
        if newFig:
            plt.figure()

        # only planet(s)' parameters
        planets_par_map = par_map[:-self.nobserv]

        P = np.array(planets_par_map[::5])
        K = np.array(planets_par_map[1::5])
        ecc = np.array(planets_par_map[2::5])
        omega = np.array(planets_par_map[3::5])
        t0 = np.array(planets_par_map[4::5])

        for planeti in range(self.nplanets):
            if self.nplanets > 1:
                # indices of the other planets
                otherplaneti = np.delete(np.arange(self.nplanets), planeti)
                # otherplaneti = int(not bool(planeti))
                otheri = otherplaneti

            t, rv, err = system.time, system.vrad, system.error # temporaries
            tt = system.get_time_to_plot(oversample=oversample)
            vel = np.zeros_like(tt)
            vel_other = np.zeros_like(t)

            # one subplot per planet
            ax = plt.subplot(self.nplanets, 1, planeti+1)

            # parameters for this planet (planeti)
            par = [P[planeti], K[planeti], ecc[planeti], omega[planeti], t0[planeti], 0.]
            # print par
            args = [tt] + par + [vel]
            get_rvn(*args)
            phase = ((tt - t0[planeti]) / P[planeti]) % 1.0

            # plot the MAP curve for this planet only (planeti)
            ax.plot(np.sort(phase), vel[np.argsort(phase)], '-g', lw=2.5, label='MAP')
            ax.plot(np.sort(phase)+1, vel[np.argsort(phase)], '-g', lw=2.5)
            ax.plot(np.sort(phase)-1, vel[np.argsort(phase)], '-g', lw=2.5)

            # plot GP predictions
            # if self.gp and plot_gp:
            #     phase = ((self.pred_t - t0) / P) % 1.0
            #     indices = np.argsort(phase)
            #     plt.plot(np.sort(phase), self.pred_y[indices], '-k', lw=0.5, alpha=0.6, label='GP mean')
                # plt.fill_between(np.sort(phase), 
                #                  y1=self.pred_y[indices]-2*self.pred_std[indices], 
                #                  y2=self.pred_y[indices]+2*self.pred_std[indices],
                #                  color='k', alpha=0.3, label='2*std')

            # the curves for the other planets (otherplaneti)
            if self.nplanets > 1:
                par = [P[otheri], K[otheri], ecc[otheri], omega[otheri], t0[otheri], 0.]
                args = [t] + par + [vel_other]
                get_rvn(*args)
            else:
                vel_other = np.zeros_like(t)

            # plot each files' values
            for i, (fname, [n, nout]) in enumerate(sorted(system.provenance.iteritems())):
                m = n-nout # how many values are there after restriction

                phase = ((t[:m] - t0[planeti]) / P[planeti]) % 1.0
                ax.errorbar(np.sort(phase), rv[np.argsort(phase)] - vsys[i] - vel_other[np.argsort(phase)],
                             yerr=err[np.argsort(phase)],
                             fmt='o'+colors[i], label=os.path.basename(fname))
                ax.errorbar(np.sort(phase)+1, rv[np.argsort(phase)] - vsys[i] - vel_other[np.argsort(phase)],
                             yerr=err[np.argsort(phase)],
                             fmt='o'+colors[i])
                ax.errorbar(np.sort(phase)-1, rv[np.argsort(phase)] - vsys[i] - vel_other[np.argsort(phase)],
                             yerr=err[np.argsort(phase)],
                             fmt='o'+colors[i])
                t, rv, err = t[m:], rv[m:], err[m:]

            # plot systematic velocity
            # for i, v in enumerate(vsys):
            #     ax.axhline(y=v, ls='--', color=colors[i], alpha=0.3)

            ax.set_xlim([-0.2, 1.2])
            ax.set_xlabel('Phase (P=%5.2f)' % P[planeti])
            ax.set_ylabel('RV [%s]'%system.units)

        if legend: plt.legend()
        plt.tight_layout()
        plt.ticklabel_format(useOffset=False)
        
        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)


    def do_plot_map_pretty(self, system, legend=True, save=None, oversample=10):

        figwidth = 3.543311946  # in inches = \hsize = 256.0748pt
        figheight = 0.75 * figwidth
        # colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        # this is Seaborn's "colorblind" pallete
        # colors = ['#0072b2', '#009e73', '#d55e00', '#cc79a7', '#f0e442', '#56b4e9']
        # this is Seaborn's "muted" pallete
        colors = ['#4878cf', '#6acc65', '#d65f5f', '#b47cc7', '#c4ad66', '#77bedb']

        t, rv, err = system.time, system.vrad, system.error # temporaries
        # tt = system.get_time_to_plot(P=self.par_map[0])
        tt = system.get_time_to_plot(oversample=oversample)
        vel = np.zeros_like(tt)
        velt = np.zeros_like(t)

        # if only one observatory, observ is always 1
        observ = np.ones_like(velt, dtype=int)
        # else it takes a different integer value for each observatory
        if self.nobserv > 1:
            n_each_observ = [v[0] for v in system.provenance.values()]
            chunks = np.split(observ, np.cumsum(n_each_observ))[:-1]
            for i in range(self.nobserv):
                chunks[i] = chunks[i] * (i+1)

            observ = np.concatenate(chunks)

        # get the RV offsets
        j = 5 if self.gp else 0
        if self.nobserv > 1:
            vsys = self.par_map[-self.nobserv:]
        else:
            vsys = [self.par_map[-1-j]]

        # we add the velocity offsets here
        # residuals will be RV-velt
        for i in range(self.nobserv):
            # ind = np.where(observ==(i+1))[0]
            # vel[ind] += vsys[i]
            ind = np.where(observ==(i+1))[0]
            velt[ind] += vsys[i]

        ## MAP estimate of the parameters
        if self.gp:
            par_map = self.par_map[:-5]
            hyper_map = self.par_map[-5:]
            # print par_map, hyper_map
            pred = gp_predictor(t, rv, err, par_map, hyper_map, 'stub')
        else:
            par_map = self.par_map

        with plt.rc_context(params):

            class MyFormatter(ticker.ScalarFormatter):
                def __call__(self, x, pos=None):
                    # call the original LogFormatter
                    rv = ticker.ScalarFormatter.__call__(self, x, pos)
                    # check if we really use TeX
                    if plt.rcParams["text.usetex"]:
                        # if we have the string ^{- there is a negative exponent
                        # where the minus sign is replaced by the short hyphen
                        rv = re.sub(r'-', r'\mhyphen', rv)
                    return rv


            fig = plt.figure(figsize=(figwidth, figheight))

            gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1], sharex=ax1)

            if self.gp and not self.only_vsys:
                # plot best solution
                # only planet(s)' parameters
                planets_par_map = self.par_map[:-5-self.nobserv]

                P = planets_par_map[::5]
                K = planets_par_map[1::5]
                ecc = planets_par_map[2::5]
                omega = planets_par_map[3::5]
                t0 = planets_par_map[4::5]
                par = [P, K, ecc, omega, t0, 0.]

                args = [tt] + par + [vel]
                get_rvn(*args)
                ax1.plot(tt, vel+vsys, '-g', lw=2.5, label='MAP')            
            elif self.only_vsys:
                pass
            else:
                # plot best solution
                # only planet(s)' parameters
                planets_par_map = par_map[:-self.nobserv]

                P = planets_par_map[::5]
                K = planets_par_map[1::5]
                ecc = planets_par_map[2::5]
                omega = planets_par_map[3::5]
                t0 = planets_par_map[4::5]
                par = [P, K, ecc, omega, t0, 0.]

                args = [tt] + par + [vel]
                get_rvn(*args)
                ax1.plot(tt, vel+vsys[0], '-k', lw=1, label='MAP')
                # for t in t0:
                # t0s = t0 + P * np.arange(-150, 51)
                # tts = get_tt(P[0], ecc[0], omega[0], t0s)
                # ax1.vlines(x=t0s, ymin=-500, ymax=500, color='k', alpha=0.8)
                # ax1.vlines(x=tts, ymin=-500, ymax=500, color='r', alpha=0.8)

            # plot GP predictions
            if self.gp:
                # ax1.plot(t, pred, 'ro')
                pred_y = self.pred_y+vsys[0]
                ax1.plot(self.pred_t, pred_y, '-k', lw=1.5, label='GP mean')
                ax1.fill_between(self.pred_t, y1=pred_y-2*self.pred_std, y2=pred_y+2*self.pred_std,
                                              color='k', alpha=0.3, label='2*std')

            # vel = np.zeros_like(t)
            if self.gp or self.only_vsys:
                # vel = self.par_map[0] * np.ones_like(tt)
                pass
            else:
                args = [t] + par + [velt]
                get_rvn(*args)

            for i, (fname, [n, nout]) in enumerate(sorted(system.provenance.iteritems())):
                m = n-nout # how many values are there after restriction
                
                # plot each files' values offset by systematic velocities
                if self.gp:
                    # here we don't remove vsys yet but this needs to be taken care of!!
                    ax1.errorbar(t[:m], rv[:m], yerr=err[:m], 
                                 fmt='o', color=colors[i], 
                                 mec='none', ms=2, capsize=0, elinewidth=0.5,
                                 label=os.path.basename(fname))
                else:
                    # because we add each offset to the model RV, we subtract them from the observations
                    ax1.errorbar(t[:m], rv[:m], yerr=err[:m], 
                                 fmt='o', color=colors[i], 
                                 mec='none', ms=2, capsize=0, elinewidth=0.5,
                                 label=os.path.basename(fname))
                # plot residuals
                if self.gp:
                    ax2.errorbar(t[:m], rv[:m]-pred[:m], yerr=err[:m], 
                                 fmt='o', color=colors[i], 
                                 mec='none', ms=2, capsize=0, elinewidth=0.5,
                                 label=fname)
                    pred = pred[m:]
                else:
                    ax2.errorbar(t[:m], rv[:m]-velt[:m], yerr=err[:m], 
                                 fmt='o', color=colors[i], 
                                 mec='none', ms=2, capsize=0, elinewidth=0.5,
                                 label=fname)
                    velt = velt[m:]

                t, rv, err = t[m:], rv[m:], err[m:]



            # plot systematic velocity
            for i, v in enumerate(vsys):
                ax1.axhline(y=v, ls='--', color=colors[i], alpha=0.6)
            ax2.axhline(y=0., ls='--', color='k', alpha=0.6)

            ax2.set_xlabel('Time [days]')
            ax1.set_ylabel('RV [%s]'%system.units)
            ax2.set_ylabel('Residuals [%s]'%system.units)

            ax1.margins(xmargin=0.05)
            ax1.minorticks_on()
            ax1.yaxis.set_major_formatter(MyFormatter())
            ax2.yaxis.set_major_formatter(MyFormatter())
            if legend: ax1.legend()
            fig.tight_layout()
            # fig.ticklabel_format(useOffset=False)

            if save:
                msg = yellow('INFO: ') + 'Saving figure to %s' % save
                clogger.info(msg)
                fig.savefig(save)


    def do_plot_map_phased_pretty(self, system, legend=False, plot_gp=True, save=None, oversample=10, title=None):
        """ 
        Phased plots ready for publication. 
        legend: put a legend in the plot, default False
        plot_gp: 
        save: save filename with extension, default None
        oversample: oversample the RV curve, default 10
        title: list of titles for each subplot (each planet), default try to guess from system's filename
        """
        # if systematic velocity only, there is nothing to do here
        if self.only_vsys: 
            msg = yellow('INFO: ') + 'This model only has offsets, a phased plot does not apply.'
            clogger.info(msg)
            return

        from math import floor, ceil
        star = get_star_name(system)

        figwidth = 3.543311946  # in inches = \hsize = 256.0748pt
        figheight = 0.75 * figwidth
        # this is Seaborn's "colorblind" pallete
        # colors = ['#0072b2', '#009e73', '#d55e00', '#cc79a7', '#f0e442', '#56b4e9']
        # this is Seaborn's "muted" pallete
        colors = ['#4878cf', '#6acc65', '#d65f5f', '#b47cc7', '#c4ad66', '#77bedb']

        # colors = 'bgrcmykw' # lets hope for less than 9 data-sets
        t, rv, err = system.time, system.vrad, system.error # temporaries
        tt = system.get_time_to_plot(oversample=oversample)
        vel = np.zeros_like(tt)

        # get the RV offsets
        j = 5 if self.gp else 0
        if self.nobserv > 1:
            vsys = self.par_map[-self.nobserv:]
        else:
            vsys = [self.par_map[-1-j]]

        ## MAP estimate of the parameters
        if self.gp:
            par_map = self.par_map[:-4]  # don't care about the hyperparameters for now
        else:
            par_map = self.par_map

        with plt.rc_context(params):

            class MyFormatter(ticker.ScalarFormatter):
                def __call__(self, x, pos=None):
                    # call the original LogFormatter
                    rv = ticker.ScalarFormatter.__call__(self, x, pos)
                    # check if we really use TeX
                    if plt.rcParams["text.usetex"]:
                        # if we have the string ^{- there is a negative exponent
                        # where the minus sign is replaced by the short hyphen
                        rv = re.sub(r'-', r'\mhyphen', rv)
                    return rv


            plt.figure(figsize=(figwidth, figheight))

            # only planet(s)' parameters
            planets_par_map = par_map[:-self.nobserv]

            P = np.array(planets_par_map[::5])
            K = np.array(planets_par_map[1::5])
            ecc = np.array(planets_par_map[2::5])
            omega = np.array(planets_par_map[3::5])
            t0 = np.array(planets_par_map[4::5])

            for planeti in range(self.nplanets):
                if self.nplanets > 1:
                    # indices of the other planets
                    otherplaneti = np.delete(np.arange(self.nplanets), planeti)
                    # otherplaneti = int(not bool(planeti))
                    otheri = otherplaneti

                t, rv, err = system.time, system.vrad, system.error # temporaries
                tt = system.get_time_to_plot(oversample=oversample)
                vel = np.zeros_like(tt)
                vel_other = np.zeros_like(t)

                # one subplot per planet
                ax = plt.subplot(self.nplanets, 1, planeti+1)
                if title is None:
                    ax.set_title(star + ' (%s)' % 'abcdefg'[planeti], loc='right', fontsize=params['font.size'])
                else:
                    ax.set_title(title[planeti], loc='right', fontsize=params['font.size'])

                # parameters for this planet (planeti)
                par = [P[planeti], K[planeti], ecc[planeti], omega[planeti], t0[planeti], 0.]
                # print par
                args = [tt] + par + [vel]
                get_rvn(*args)
                phase = ((tt - t0[planeti]) / P[planeti]) % 1.0

                # plot the MAP curve for this planet only (planeti)
                ax.plot(np.sort(phase), vel[np.argsort(phase)], '-', color='k', lw=1, label='MAP')
                ax.plot(np.sort(phase)+1, vel[np.argsort(phase)], '-', color='k', lw=1)
                ax.plot(np.sort(phase)-1, vel[np.argsort(phase)], '-', color='k', lw=1)

                # plot GP predictions
                # if self.gp and plot_gp:
                #     phase = ((self.pred_t - t0) / P) % 1.0
                #     indices = np.argsort(phase)
                #     plt.plot(np.sort(phase), self.pred_y[indices], '-k', lw=0.5, alpha=0.6, label='GP mean')
                    # plt.fill_between(np.sort(phase), 
                    #                  y1=self.pred_y[indices]-2*self.pred_std[indices], 
                    #                  y2=self.pred_y[indices]+2*self.pred_std[indices],
                    #                  color='k', alpha=0.3, label='2*std')

                # the curves for the other planets (otherplaneti)
                if self.nplanets > 1:
                    par = [P[otheri], K[otheri], ecc[otheri], omega[otheri], t0[otheri], 0.]
                    args = [t] + par + [vel_other]
                    get_rvn(*args)
                else:
                    vel_other = np.zeros_like(t)

                # plot each files' values
                for i, (fname, [n, nout]) in enumerate(sorted(system.provenance.iteritems())):
                    m = n-nout # how many values are there after restriction

                    phase = ((t[:m] - t0[planeti]) / P[planeti]) % 1.0
                    points = rv[np.argsort(phase)] - vsys[i] - vel_other[np.argsort(phase)]
                    e = err[np.argsort(phase)]

                    # limits for plot
                    max_value, min_value = max(points+e), min(points-e)
                    max_y = ceil(max(abs(max_value), abs(min_value)))
                    print max_y
                    min_y = -max_y

                    ax.errorbar(np.sort(phase), points,
                                 yerr=err[np.argsort(phase)],
                                 fmt='o', color=colors[i], 
                                 mec='none', ms=2, capsize=0, elinewidth=0.5,
                                 label=os.path.basename(fname))
                    ax.errorbar(np.sort(phase)+1, rv[np.argsort(phase)] - vsys[i] - vel_other[np.argsort(phase)],
                                 yerr=err[np.argsort(phase)],
                                 fmt='o', color=colors[i], 
                                 mec='none', ms=2, capsize=0, elinewidth=0.5)
                    ax.errorbar(np.sort(phase)-1, rv[np.argsort(phase)] - vsys[i] - vel_other[np.argsort(phase)],
                                 yerr=err[np.argsort(phase)],
                                 fmt='o', color=colors[i], 
                                 mec='none', ms=2, capsize=0, elinewidth=0.5)
                    t, rv, err = t[m:], rv[m:], err[m:]

                # plot systematic velocity
                # for i, v in enumerate(vsys):
                #     ax.axhline(y=v, ls='--', color=colors[i], alpha=0.3)

                ax.set_xlim([-0.1, 1.1])
                # ax.set_ylim([min_y, max_y])
                ax.margins(ymargin=0.05)
                ax.set_xlabel(r'$\phi$ (P = %3.2f days)' % P[planeti])
                ax.set_ylabel('RV [%s]'%system.units)
                ax.minorticks_on()
                ax.xaxis.set_major_formatter(MyFormatter())
                ax.yaxis.set_major_formatter(MyFormatter())

            if legend: plt.legend()
            plt.tight_layout()
            
            if save:
                msg = yellow('INFO: ') + 'Saving figure to %s' % save
                clogger.info(msg)
                plt.savefig(save)


    def do_plot_map_and_ml(self, system, legend=True, save=None):
        colors = 'kbgrcmyw' # lets hope for less than 9 data-sets
        t, rv, err = system.time, system.vrad, system.error # temporaries
        # tt = system.get_time_to_plot(P=self.par_map[0])
        tt = system.get_time_to_plot()
        vel_map = np.zeros_like(tt)
        vel_mle = np.zeros_like(tt)

        nobserv = len(system.provenance)  # number of observatories
        j = 4 if self.gp else 0
        if nobserv > 1:
            map_vsys = self.par_map[-nobserv+1:]
            mle_vsys = self.par_mle[-nobserv+1:]
            get_rvn = get_rvn_ms
        else:
            map_vsys = self.par_map[-1-j]
            mle_vsys = self.par_mle[-1-j]
            get_rvn = get_rvn_os

        ## MAP estimate of the parameters
        if self.only_vsys:  # systematic velocity only
            vel_map = self.par_map[0] * np.ones_like(tt)
        elif self.gp:
            par_map = self.par_map[:-4] 
            hyper_map = self.par_map[-4:]
            if self.gp_only:
                pred = gp_predictor(t, rv, err, par_map, hyper_map, 'constant')
            else:
                pred = gp_predictor(t, rv, err, par_map, hyper_map, 'keplerian')
        else:
            par_map = self.par_map

        ## ML estimate of the parameters
        if self.only_vsys:  # systematic velocity only
            vel_mle = self.par_mle[0] * np.ones_like(tt)
        elif self.gp:
            par_mle = self.par_mle[:-4] 
            hyper_mle = self.par_mle[-4:]
            if self.gp_only:
                pred = gp_predictor(t, rv, err, par_mle, hyper_mle, 'constant')
            else:
                pred = gp_predictor(t, rv, err, par_mle, hyper_mle, 'keplerian')
        else:
            par_mle = self.par_mle


        newFig=True
        if newFig:
            fig = plt.figure()

        gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharex=ax1)

        if self.only_vsys:  # systematic velocity only
            ax1.plot(tt, vel_map, '-g', lw=2.5, label='MAP')
            ax1.plot(tt, vel_mle, '--r', lw=2.5, label='ML')
        elif not self.gp_only:
            # plot best solution
            map_P = par_map[:-1:5]
            map_K = par_map[1:-1:5]
            map_ecc = par_map[2:-1:5]
            map_omega = par_map[3:-1:5]
            map_t0 = par_map[4:-1:5]
            par = [map_P, map_K, map_ecc, map_omega, map_t0, map_vsys]

            if nobserv > 1:
                observ = np.ones_like(vel_map)
                args = [tt] + par + [observ] + [vel_map]
            else:
                args = [tt] + par + [vel_map]
            get_rvn(*args)
            ax1.plot(tt, vel_map, '-g', lw=2.5, label='MAP')

            mle_P = par_mle[:-1:5]
            mle_K = par_mle[1:-1:5]
            mle_ecc = par_mle[2:-1:5]
            mle_omega = par_mle[3:-1:5]
            mle_t0 = par_mle[4:-1:5]
            par = [mle_P, mle_K, mle_ecc, mle_omega, mle_t0, mle_vsys]

            if nobserv > 1:
                observ = np.ones_like(vel_mle)
                args = [tt] + par + [observ] + [vel_mle]
            else:
                args = [tt] + par + [vel_mle]
            get_rvn(*args)
            ax1.plot(tt, vel_mle, '--r', lw=2.5, label='ML')



        # plot GP predictions
        if self.gp:
            ax1.plot(self.pred_t, self.pred_y, '-k', lw=1.5, label='GP mean')
            ax1.fill_between(self.pred_t, y1=self.pred_y-2*self.pred_std, y2=self.pred_y+2*self.pred_std,
                                          color='k', alpha=0.3, label='2*std')

        vel = np.zeros_like(t)
        if self.gp_only or self.only_vsys:
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


        # plot systematic velocity
        ax1.axhline(y=map_vsys, ls='--', color='k', alpha=0.3)

        ax2.set_xlabel('Time [days]')
        ax1.set_ylabel('RV [%s]'%system.units)
        ax2.set_ylabel('Residuals [%s]'%system.units)
        if legend: ax1.legend()
        fig.tight_layout()
        # fig.ticklabel_format(useOffset=False)

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            fig.savefig(save)


    def do_hist_plots(self, show_priors=True, save=None):
        try:
            self.trace
        except AttributeError:
            self.read_iterations()       

        j = 1 if self.jitter else 0
        j = j+4 if self.gp else j+0

        if show_priors: 
            from .prior_funcs import uniform, jeffreys, modjeffreys, beta
            # print 
            # msg = blue('INFO: ') + 'Plotting logP. prior(P)=Jeffreys, prior(logP)=Uniform'
            # clogger.info(msg)

        # Sturges' formula for the optimal number of bins in the histograms
        n = self.nest_samples.shape[1]
        k = int(np.log2(n) + 1)
        k = 100

        fig_ysize = 6 if self.npar==2 else 12
        fig_ysize = 4 if self.only_vsys else fig_ysize


        plt.figure()
        nsubplots = 1
        ax1 = plt.subplot(j+1,1,nsubplots)
        nsubplots += 1
        # ax1.hist(self.nest_samples[-1, burnin:].T, bins=k, normed=True, label='vsys')
        # hist(self.ptrace.par7, bins='knuth', ax=ax1, normed=True, color='r')
        hist(self.nest_samples[-j, burnin:].T, bins='knuth', ax=ax1, normed=True, label='vsys', alpha=0.5)
        ax1.set_xlabel('Vsys [m/s]')

        if self.jitter:
            ax2 = plt.subplot(j+1,1,nsubplots)
            nsubplots += 1
            # ax2.hist(self.nest_samples[-2, burnin:].T, bins=k, normed=True)
            # hist(self.ptrace.par6, bins='knuth', ax=ax2, normed=True, color='r')
            hist(self.nest_samples[-j-1, burnin:].T, bins='knuth', ax=ax2, normed=True, alpha=0.5)
            ax2.set_xlabel('jitter [m/s]')

        if self.trend:
            ax2 = plt.subplot(j+1,1,nsubplots)
            nsubplots += 1

            # tr = np.loadtxt(self.root+'.txt', usecols=(3,))
            # hist(tr, color='r', bins=100, ax=ax2, normed=True)

            # tr = np.loadtxt(self.root+'post_equal_weights.dat', usecols=(1,))
            # hist(tr, color='g', bins=100, ax=ax2, normed=True)            

            # ax2.hist(self.nest_samples[-2, burnin:].T, bins=k, normed=True)
            hist(self.nest_samples[j, burnin:].T, bins=100, ax=ax2, normed=True)
            # ax2.axvline(x=self.par_map[j], lw=2, color='k')
            # ax2.axvline(x=self.par_mle[j], lw=2, color='k')
            # ax2.axvline(x=self.par_mean[j], lw=2, color='k')
            ax2.set_xlabel('trend [m/s]')


        plt.tight_layout()
        if self.nplanets == 0: return

        plt.figure(figsize=(8, fig_ysize))
        gs = gridspec.GridSpec(min(5,self.npar), 1)

        # period histogram(s)
        ax1 = plt.subplot(gs[0])
        # ax1.hist(self.nest_samples[0:-1-j:5, burnin:].T.flatten(), bins=100, range=[50, 80], normed=True, label=['planet1', 'planet2', 'planet3'])
        hist(self.ptrace.par1, ax=ax1, bins='knuth', normed=True, color='r')
        hist(self.nest_samples[0:-1-j:5, :].T.flatten(), ax=ax1, bins='knuth', normed=True, alpha=0.5, label=['planet1', 'planet2', 'planet3'])
        # if show_priors:
        #     xx = np.logspace(np.log10(1.), np.log10(4000), 200)
        #     ax1.plot(xx, jeffreys(xx, 1., 4000), lw=3, color='k')
        ax1.set_xlabel('P [d]')
        ax1.legend(frameon=False)

        # semi amplitude histogram(s) - actually log of it because usually the prior spans a big range
        ax2 = plt.subplot(gs[1])
        ax2.hist(self.nest_samples[1:-1-j:5, burnin:].T, bins=200, normed=True, range=[0, 10])
        if show_priors:
        #     m = np.min(self.nest_samples[1:-1-j:5, :].T)
        #     xx = np.logspace(np.log(m), np.log10(2129), 200)
        #     ax2.plot(np.log(xx), modjeffreys(xx, 1., 2129), lw=3, color='k')
            xx = np.linspace(0, 60, 200)
            ax2.plot(xx, modjeffreys(xx, 1., 60), lw=3, color='k')
        ax2.set_xlabel('K [m/s]')

        # eccentricity histogram(s)
        ax3 = plt.subplot(gs[2])
        ax3.hist(self.posterior_samples[2:-1-j:5, :].T, bins=k, normed=True)
        if show_priors:
            # xx = np.linspace(0, 1, 200)
            # ax3.plot(xx, uniform(xx, 0, 1), lw=3, color='k')
            xx = np.linspace(0.0001, 1, 500)
            ax3.plot(xx, beta(xx, 0.867, 3.03), lw=3, color='k')
        ax3.set_xlabel('eccentricity')

        # omega histogram(s)
        ax4 = plt.subplot(gs[3])
        ax4.hist(self.posterior_samples[3:-1-j:5, :].T, bins=k, normed=True)
        if show_priors:
            xx = np.linspace(0, 2*np.pi, 200)
            ax4.plot(xx, uniform(xx, 0, 2*np.pi), lw=3, color='k')
        ax4.set_xlim([0., 2*np.pi])
        ax4.set_xlabel(r'$\omega$')

        # t0 histogram(s)
        ax5 = plt.subplot(gs[4])
        ax5.hist(self.posterior_samples[4:-1-j:5, :].T, bins=k, normed=True)
        ax5.set_xlabel(r'$t_0$')


        plt.tight_layout()

        if save:
            msg = yellow('INFO: ') + 'Saving figure to %s' % save
            clogger.info(msg)
            plt.savefig(save)

    def do_triangle_plot(self, planet='both'):
        """
        Triangle plot of posterior samples
        """
        try:
            self.trace
        except AttributeError:
            self.read_iterations()

        triangle_plot(self.posterior_samples[:5,:].T, quantiles=[0.5])

    def do_triangle_plot_kde(self):

        triangle_plot_kde(self.posterior_samples[:5,:].T)



