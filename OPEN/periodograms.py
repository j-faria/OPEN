# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
from classes import PeriodogramBase
from numpy import *
import numpy as np
from numpy.fft import *  
from scipy.stats import rankdata
import matplotlib.pyplot as plt
from sys import float_info
from multiprocessing import cpu_count

from .shell_colors import yellow, blue, red
from .logger import clogger, logging
from .utils import ask_yes_no, get_star_name

from ext.blombscargle import blombscargle
from ext.glombscargle import glombscargle, glombscargle_extra_out


try:
  from gatspy import periodic
  gatspy_is_available = True
except ImportError:
  gatspy_is_available = False


help_text = """
This command calculates a number of different periodograms 
of the observed radial velocities. The algorithms are:

  -g --gls
    This is the default option. 
    It calculates the Generalised Lomb-Scargle periodogram as it is
    defined in Zechmeister & Kürster, A&A 496, 577–584 (2009).
    This is equivalent to least-squares fitting of sine waves,
    includes an offset and individual weights. 
    The default normalisation is that of Horne & Baliunas, ApJ, 302, 757-763 (1986).

  -m --bgls
    This option calculates the Bayesian Generalized Lomb-Scargle (BGLS) periodogram
    as defined by Mortier et al (2014).

  -b --bayes       Calculate the Bayesian LS periodogram
    Description

  -l --ls
    Description

  -z --hoef        Calculate the Hoeffding-test "periodogram" with Zucker's algorithm
    Description

  -r --multiband   Calculate the multiband periodogram; Vanderplas & Ivezic (2015)
    Description
"""
help_text = blue(' : ').join([s+'\n' for s in help_text.split('\n')])



# This is the Lomb-Scargle periodogram calculated using the fast 
# algorithm from Press & Rybicki (1989).
#####################################################################
class ls_PressRybicki(PeriodogramBase):
  """
  Compute the Lomb-Scargle periodogram based on the algorithm by
    Press, W. H. & Rybicki, G. B. 1989 
    ApJ vol. 338, p. 277-280. 
    Fast algorithm for spectral analysis of unevenly sampled data 
    bib code: 1989ApJ...338..277P 

  Parameters
  ----------
    rv : RVSeries
        The radial velocity curve or any object providing the attributes
        time, vrad which define the data.
    ofac : int
        Oversampling factor (default=6).
    hifac : float
        hifac * "average" Nyquist frequency is highest frequency for 
        which the periodogram will be calculated (default=6).

  Attributes
  ----------
    power : array
        The periodogram power.
    freq : array
        The frequency array.
    jmax : array
        The array index corresponding to the max(power).
    fap : float
        False Alarm Probability of the largest periodogram value.
  """
  def __init__(self, rv, ofac=6, hifac=6, quantity='vrad'):
    self.name = 'Lomb-Scargle - Press & Rybicki algorithm'

    self.power = None
    self.freq = None
    self.ofac = ofac
    self.hifac = hifac
    if quantity == 'vrad':
      self.y = rv.vrad
    elif quantity == 'bis':
      self.y = rv.extras.bis_span
    elif quantity == 'fwhm':
      if ask_yes_no('Should I remove a linear trend first? (y/N) ', False):
        m, b = polyfit(self.t, rv.extras.fwhm, 1)
        yp = polyval([m, b], self.t)
        self.y = rv.extras.fwhm - yp
      else:
        self.y = rv.extras.fwhm
    elif quantity == 'rhk':
      self.y = rv.extras.rhk
    elif quantity == 'contrast':
      self.y = rv.extras.contrast
    elif quantity == 'resid':
      try:
        self.y = rv.fit['residuals']
      except TypeError:
        clogger.fatal('error!')
        return 
    self.t = rv.time
    self.th = rv.time - min(rv.time)
    self.jmax = self.fap = None

    self.__fasper()

  def __spread(self, y, yy, n, x, m):  
    """ 
    Given an array yy(0:n-1), extirpolate (spread) a value y into 
    m actual array elements that best approximate the "fictional" 
    (i.e., possible noninteger) array element number x. The weights 
    used are coefficients of the Lagrange interpolating polynomial.
    """  
    nfac=[0,1,1,2,6,24,120,720,5040,40320,362880]  
    if m > 10. :  
      print 'factorial table too small in spread'  
      return  
    
    ix=long(x)  
    if x == float(ix):   
      yy[ix]=yy[ix]+y  
    else:  
      ilo = long(x-0.5*float(m)+1.0)  
      ilo = min( max( ilo , 1 ), n-m+1 )   
      ihi = ilo+m-1  
      nden = nfac[m]  
      fac=x-ilo  
      for j in range(ilo+1,ihi+1): fac = fac*(x-j)  
      yy[ihi] = yy[ihi] + y*fac/(nden*(x-ihi))  
      for j in range(ihi-1,ilo-1,-1):  
        nden=(nden/(j+1-ilo))*(j-ihi)  
        yy[j] = yy[j] + y*fac/(nden*(x-j))  


  def __fasper(self, MACC=4):
    """
    Given abscissas self.t (which need not be equally spaced) and 
    ordinates self.y, and given a desired oversampling factor ofac,
    this routine creates an array freq with a sequence of nout 
    increasing frequencies (not angular frequencies) up to hifac 
    times the "average" Nyquist frequency, and creates an array power
    with the values of the Lomb normalized periodogram at those 
    frequencies. This routine also returns jmax such that power(jmax)
    is the maximum element in power, and prob, an estimate of the 
    significance of that maximum against the hypothesis of random 
    noise. A small value of prob indicates that a significant 
    periodic signal is present. 
    """  

    x = self.t
    y = self.y
    ofac, hifac = self.ofac, self.hifac

    #Check dimensions of input arrays  
    n = long(len(x))
    if n != len(y):  
      print 'Incompatible arrays.'  
      return  
    
    nout  = 0.5*ofac*hifac*n  
    nfreqt = long(ofac*hifac*n*MACC)  # Size the FFT as next power  
    nfreq = 64L                       # of 2 above nfreqt.  
    
    while nfreq < nfreqt:   
      nfreq = 2*nfreq  
    
    ndim = long(2*nfreq)  
      
    ## Compute the mean
    ave = y.mean()  
    ##sample variance because the divisor is N-1  
    var = ((y-y.mean())**2).sum()/(len(y)-1)   
    # and range of the data.  
    xmin = x.min()  
    xmax = x.max()  
    xdif = xmax-xmin  
    
    #extirpolate the data into the workspaces  
    wk1 = zeros(ndim, dtype='complex')  
    wk2 = zeros(ndim, dtype='complex')  
    
    fac  = ndim/(xdif*ofac)  
    fndim = ndim  
    ck  = ((x-xmin)*fac) % fndim  
    ckk  = (2.0*ck) % fndim  
    
    for j in range(0L, n):  
      self.__spread(y[j]-ave,wk1,ndim,ck[j],MACC)  
      self.__spread(1.0,wk2,ndim,ckk[j],MACC)  
    
    #Take the Fast Fourier Transforms  
    wk1 = ifft( wk1 )*len(wk1)  
    wk2 = ifft( wk2 )*len(wk1)  
    
    wk1 = wk1[1:nout+1]  
    wk2 = wk2[1:nout+1]  
    rwk1 = wk1.real  
    iwk1 = wk1.imag  
    rwk2 = wk2.real  
    iwk2 = wk2.imag  
      
    df  = 1.0/(xdif*ofac)  
      
    #Compute the Lomb value for each frequency  
    hypo2 = 2.0 * abs( wk2 )  
    hc2wt = rwk2/hypo2  
    hs2wt = iwk2/hypo2  
    
    cwt  = sqrt(0.5+hc2wt)  
    swt  = sign(hs2wt)*(sqrt(0.5-hc2wt))  
    den  = 0.5*n+hc2wt*rwk2+hs2wt*iwk2  
    cterm = (cwt*rwk1+swt*iwk1)**2./den  
    sterm = (cwt*iwk1-swt*rwk1)**2./(n-den)  
    
    wk1 = df*(arange(nout, dtype='float')+1.)  
    wk2 = (cterm+sterm)/(2.0*var)  
    pmax = wk2.max()  
    jmax = wk2.argmax()  
    
    #Estimate significance of largest peak value  
    expy = exp(-pmax)            
    effm = 2.0*(nout)/ofac         
    prob = effm*expy  
    
    if prob > 0.01:   
      prob = 1.0-(1.0-expy)**effm  
    
    #Assign local variables to class attributes
    self.freq, self.power = wk1, wk2
    self.nf = nout
    self.fap, self.jmax = prob, jmax


# This is the ...... periodogram etc etc etc
#
#
#####################################################################
class LombScargle(PeriodogramBase):
  """
    The constructor of *LombScargle* takes a *TimeSeries* instance, i.e., a 
    light curve object, as first argument. It then computes the usual 
    Lomb-Scargle periodogram using a fast algorithm. The frequency array is
    constructed on the fly based on the oversampling keywords, which are
    mandatory in this case. The power is normalized according to the
    prescription of [HB86]_.
  """
  
  def __init__(self, rv, ofac=6, hifac=1):
    """
      Parameters:
       - `rv` - TimesSeries instance, The light curve to be analyzed.
       - `ofac` - int, Oversampling factor.
       - `hifac` - float, Maximum frequency `freq` = `hifac` * (average Nyquist frequency).
      
      .. note::
        Adapted from routine of the same routine in [NR]_ , 
        based on period.pro by Han Wen, August 1996.
      
      The result, i.e., the power, is stored in the class property `power`.
    """

    self.name = 'Lomb Scargle'
    self.freq=None
    self.power=None
    self.ofac=ofac
    self.hifac=hifac
    self.t=rv.time
    self.y=rv.vrad
    self.label = {'title': 'Lomb-Scargle periodogram', \
                  'xlabel': 'Frequency', \
                  'ylabel': 'Scargle Power'}
    self.__calcPeriodogram()


  def __calcPeriodogram(self):

    ofac=int(self.ofac); hifac=float(self.hifac)

    if self.freq is None:
      self.__buildFreq(ofac,hifac)
    nout = int(ofac*hifac*len(self.t)/2.)
    vari = var(self.y)
    ave = mean(self.y)

    xdif = max(self.t)-min(self.t)
    xave = 0.5*(max(self.t)+min(self.t))

    # Get starting frequency
    pnow = 1./(xdif*ofac)

    arg = 2.*pi*(self.t-xave)*pnow
    wpr = -2.*sin(0.5*arg)**2
    wpi = sin(arg)
    wr = cos(arg)
    wi = wpi

    py = zeros(nout)
    yy = self.y-ave

    for i in range(nout):
      sumsh = sum(wi*wr)
      sumc = sum((wr-wi)*(wr+wi))
      wtau = 0.5 * arctan(2.*sumsh/sumc)
      swtau = sin(wtau)
      cwtau = cos(wtau)
      ss = wi*cwtau-wr*swtau
      cc = wr*cwtau+wi*swtau
      sums = sum(ss**2)
      sumc = sum(cc**2)
      sumsy = sum(yy*ss)
      sumcy = sum(yy*cc)
      wtemp = wr
      wr   =(wtemp*wpr-wi*wpi)+wr
      wi   =wi*wpr+wtemp*wpi+wi
      py[i]=0.5*(sumcy**2/sumc+sumsy**2/sums)/vari

    self.power = py.copy()
    self.M = 2.*nout/ofac
    self.N = len(self.y)

  def __buildFreq(self,ofac,hifac):
    nout = int(ofac*hifac*len(self.t)/2.)
    xdif = max(self.t)-min(self.t)
    pnow = 1./(xdif*ofac)
    px = pnow + arange(nout)/(ofac*xdif)
    self.freq = px
    return 0

  def prob(self, Pn):
    """
      Returns the probability to obtain a power *Pn* or larger from the noise,
      which is assumes to be Gaussian.
      
      Parameters:
        - `Pn` - float, Power threshold.

      .. note::
        *LombScargle* calculates the quantity (N-1)/2.*p=p' (in the formalism of 
        [ZK09]_), which is de facto the normalization
        prescription of [HB86]_. In this
        scheme the probability P(p'>Pn) is given by the following statement:
        
        .. math::
          P(p'>Pn) = \\left(1 - 2 \\frac{Pn}{N-1} \\right)^{(N-3)/2}

        If properly normalized to the population variance of the time series, which
        must be known a priori (usually not the case), the
        power :math:`p/p_n=p"` is a direct measure of the SNR as proposed by [Scargle82]_:
        
        .. math::
          P(p">Pn) = exp(-Pn) \\; .
        
        This formula is often used erroneously in this context.
    """
    return (1.-2.*Pn/(self.N-1.))**((self.N-3.)/2.)


# This is the Generalized Lomb-Scargle (GLS) periodogram as defined
# by Zechmeister & Kürster (2009). 
#####################################################################
class gls(PeriodogramBase):
  """
  Compute the Generalized Lomb-Scargle (GLS) periodogram.

  This class implements the error-weighted Lomb-Scargle periodogram as
  developed by [ZK09]_ using various possible normalizations.

  The constructor takes a RVSeries instance (i.e. a rv curve) as 
  first argument. As the algorithm is slow-ish, an implementation
  in Fortran is available if the keyword 'ext' is set to True in 
  the constructor or globally. 
  There is an optional `freq` array, that can contain the 
  frequencies on which to calculate the periodogram. If not provided
  (....)
    
  Parameters
    rv : RVSeries
        The radial velocity curve or any object providing the attributes
        time, vrad and error which define the data.
    ofac : int
        Oversampling factor (default=6).
    hifac : float
        hifac * "average" Nyquist frequency is highest frequency for 
        which the periodogram will be calculated (default=1).
    freq : array, optional
        Contains the frequencies at which to calculate the periodogram.
        If not given, a frequency array will be automatically generated.
    quantity : string, optional
        For which quantity to calculate the periodogram. Possibilities are
        'bis', 'rhk', 'contrast' or 'fwhm' other than the default 'vrad'.
    norm : string, optional
        The normalization; either "Scargle", "HorneBaliunas", or 
        "Cumming". Default is "HorneBaliunas".
    stats : boolean, optional
        Set True to obtain some statistical output (default is False).
    ext : boolean, optional
        Use Fortran extension in the calculation (default is True)
  
  Attributes
    power : array
        The normalized power of the GLS.
    freq : array
        The frequency array.
    ofac : int
        The oversampling factor.
    hifac : float
        The maximum frequency.
    norm : string
        The normalization used.
  """

  def __init__(self, rv, ofac=6, hifac=1, freq=None, quantity='vrad',
               norm="HorneBaliunas", stats=False, ext=True, force_notrend=False, full_output=False):
    self.name = 'Generalized Lomb-Scargle'
    try:
      self.star_name = get_star_name(rv)
    except AttributeError:
      self.star_name = ''

    self.power = None
    self.freq = freq
    self.ofac, self.hifac = ofac, hifac
    self.fullout = full_output
    self.t = rv.time
    self.th = rv.time - min(rv.time)
    if quantity == 'vrad':
      self.y = rv.vrad
      self.error = rv.error
    elif quantity == 'bis':
      self.y = rv.extras.bis_span
      self.error = 2. * rv.error #ones_like(self.y)
    elif quantity == 'fwhm':
      if not force_notrend and ask_yes_no('Should I remove a linear trend first? (y/N) ', False):
        m, b = polyfit(self.t, rv.extras.fwhm, 1)
        yp = polyval([m, b], self.t)
        self.y = rv.extras.fwhm - yp
      else:
        self.y = rv.extras.fwhm
      self.error = 2.35 * rv.error #ones_like(self.y)
    elif quantity == 'rhk':
      self.y = rv.extras.rhk
      self.error = rv.extras.sig_rhk
    elif quantity == 'contrast':
      self.y = rv.extras.contrast
      self.error = ones_like(self.y)
    elif quantity == 'resid':
      try:
        self.y = rv.fit['residuals']
      except TypeError:
        clogger.fatal('error!')
        return 
      self.error = rv.error

    self.norm = norm
    # Check and assign normalization
    self.label = {'title': 'Generalized Lomb Periodogram',\
                  'xlabel': 'Frequency',\
                  'ylabel': "Normalization not implemented!"}
    if self.norm == "Scargle":
      self.label["ylabel"] = "Normalized Power (Scargle 1982)"
    elif self.norm == "HorneBaliunas":
      self.label["ylabel"] = "Normalized Power (Horne & Baliunas 1986)"
    elif self.norm == "Cumming":
      self.label["ylabel"] = "Normalized Power (Cumming 1999)"
    else:
      pass
    
    self._stats = stats
    if ext: 
      self.__calcPeriodogramFast()
    else: 
      self.__calcPeriodogram()
    
  def __calcPeriodogram(self):
    """ Compute the GLS. Notation very close to [ZK09]_ """

    # Build frequency array if not present
    if self.freq is None:
      plow = max(0.5, 2*ediff1d(self.t).min()) # minimum meaningful period?
      self.__buildFreq(plow=plow)
    # Circular frequencies
    omegas = 2.*pi * self.freq

    err2 = self.error*self.error
    w = (1. / sum(1./err2)) / err2

    self._Y = sum(w*self.y)        # Eq. (7)
    yh = self.y - self._Y          # Subtract weighted mean
    self._YY = sum(w * yh**2)
    
    # Unnormalized power
    self._upow = zeros(len(omegas))
    self._a = zeros(len(omegas))
    self._b = zeros(len(omegas))
    self._off = zeros(len(omegas))
    for i, omega in enumerate(omegas):
      x = omega*self.th
      cosx = cos(x)
      sinx = sin(x)
      wcosx = w*cosx         # attach weights
      wsinx = w*sinx         # attach weights
      
      C = sum(wcosx)         # Eq. (8)
      S = sum(wsinx)         # Eq. (9)

      YC = sum(yh*wcosx)     # Eq. (11)
      YS = sum(yh*wsinx)     # Eq. (12)
      CCh = sum(wcosx*cosx)  # Eq. (13)
      CSh = sum(wcosx*sinx)  # Eq. (15)
      SSh = 1.-CCh
      CC = CCh-C*C           # Eq. (13)
      SS = SSh-S*S           # Eq. (14)
      CS = CSh-C*S           # Eq. (15)
      D = CC*SS-CS*CS        # Eq. (6)
      
      self._a[i] = (YC*SS-YS*CS) / D
      self._b[i] = (YS*CC-YC*CS) / D
      self._off[i] = -self._a[i]*C - self._b[i]*S

      self._upow[i] = (SS*YC*YC + CC*YS*YS - 2.*CS*YC*YS) / (self._YY*D) # Eq. (5) in ZK09

    self.N = len(self.y)
    # An ad-hoc estimate of the number of independent frequencies 
    # see discussion following Eq. (24)
    self.M = (max(self.freq)-min(self.freq)) * (self.th.max() - self.th.min())

    # Normalization:
    if self.norm == "Scargle":
      popvar=raw_input('Input a priori known population variance:')
      self.power = self._upow/float(popvar)
    if self.norm == "HorneBaliunas":
      self.power = (self.N-1.)/2.*self._upow
    if self.norm == "Cumming":
      self.power = (self.N-3.)/2. * self._upow/(1.-max(self._upow))
    
    # Output statistics
    if self._stats:
      self._output()
    # if self._showPlot:
    #   self._plot()

  def __calcPeriodogramFast(self):
    """ Compute the GLS using the Fortran extension 
    which allows speedups of ~6x."""
    ncpu = cpu_count()

    # Build frequency array if not present
    if self.freq is None:
      plow = max(0.5, 2*ediff1d(self.t).min()) # minimum meaningful period?
      self.plow = plow
      self.__buildFreq(plow=plow)
    # Circular frequencies
    omegas = 2.*pi * self.freq
    self.omegas = omegas

    if self.fullout:
      self._upow, self.M, self.A, self.B, self.offset = glombscargle_extra_out(self.t, self.y, self.error, omegas, ncpu)
    else:
      # unnormalized power and an estimate of the number of independent frequencies 
      self._upow, self.M = glombscargle(self.t, self.y, self.error, omegas, ncpu)

    self.N = len(self.y)
    # Normalization:
    self._normalize()
    
    # Output statistics
    if self._stats:
      self._output()
    # if self._showPlot:
    #   self._plot()

  def _normalize(self):
    if self.norm == "Scargle":
      popvar=raw_input('pyTiming::gls - Input a priori known population variance:')
      self.power = self._upow/float(popvar)
    if self.norm == "HorneBaliunas":
      self.power = (self.N-1.)/2.*self._upow
    if self.norm == "Cumming":
      self.power = (self.N-3.)/2. * self._upow/(1.-max(self._upow))

  def _normalize_value(self, a):
    # normalize the value (or array) a
    if self.norm == "Scargle":
      popvar=raw_input('pyTiming::gls - Input a priori known population variance:')
      return a/float(popvar)
    if self.norm == "HorneBaliunas":
      return (self.N-1.)/2.*a
    if self.norm == "Cumming":
      return (self.N-3.)/2. * a/(1.-max(np.atleast_1d(a)))    

  def __buildFreq(self, plow=0.5):
    """
      Build frequency array (`freq` attribute).
    """
    xdif = max(self.th)-min(self.th)
    # print xdif
    # nout = self.ofac * self.hifac * len(self.th)/2
    nout = int(xdif * self.ofac / plow)
    # print self.ofac, nout
    # sys.exit(0)
    self.freq = 1./(xdif) + arange(nout)/(self.ofac*xdif)

  def prob(self, Pn):
    """
      Probability of obtaining the given power.
    
      Calculate the probability to obtain a power higher than
      `Pn` from the noise, which is assumed to be Gaussian.
      
      .. note:: This depends on the normalization
        (see [ZK09]_ for further details).

        - `Scargle`: 
        .. math::
          exp(-Pn)

        - `HorneBaliunas`: 
        .. math::
          \\left(1 - 2 \\times \\frac{Pn}{N-1} \\right)^{(N-3)/2}
        
        - `Cumming`: 
        .. math:: 
          \\left(1+2\\times \\frac{Pn}{N-3}\\right)^{-(N-3)/2}
      
      Parameters
        Pn : float
          Power threshold.
      
      Returns
        Probability : float
          The probability to obtain a power equal or
          higher than the threshold from the noise.
    """
    if self.norm=="Scargle": return exp(-Pn)
    if self.norm=="HorneBaliunas": return (1.-2.*Pn/(self.N-1.))**((self.N-3.)/2.)
    if self.norm=="Cumming": return (1.+2.*Pn/(self.N-3.))**(-(self.N-3.)/2.)

  def probInv(self, Prob):
    """
      Calculate minimum power for a given probability.
    
      This function is the inverse of `Prob(Pn)`.
      Returns the minimum power for a given probability threshold Prob.
      
      Parameters
        Prob : float
          Probability threshold.
      
      Returns
        Power threshold : float
          The minimum power for the given
          false-alarm probability threshold.
    """
    if self.norm=="Scargle": return -log(Prob)
    if self.norm=="HorneBaliunas": return (self.N-1.)/2.*(1.-Prob**(2./(self.N-3.)))
    if self.norm=="Cumming": return (self.N-3.)/2.*(Prob**(-2./(self.N-3.))-1.)



# This is the Multiband Lomb-Scargle (GLS) periodogram as introduced
# by Vanderplas & Ivezic, ApJ 812:18 (2015)
#####################################################################
class MultiBandGLS(PeriodogramBase):
  """
  Compute the Multiband Lomb-Scargle periodogram.
  We use this periodogram when there is data after the HARPS fiber update.
  This depends on gatspy (http://www.astroml.org/gatspy)

  The constructor takes a RVSeries instance (i.e. a rv curve) as first argument.
  By default the time of offset is May 28, 2015 the date of the HARPS upgrade,
  but this can be changed by providing the time_of_offset argument.
    
  Parameters
    rv : RVSeries
        The radial velocity curve or any object providing the attributes
        time, vrad and error which define the data.
    ofac : int
        Oversampling factor (default=6).
    hifac : float
        hifac * "average" Nyquist frequency is highest frequency for 
        which the periodogram will be calculated (default=1).
    quantity : string, optional
        For which quantity to calculate the periodogram. Possibilities are
        'bis', 'rhk', 'contrast' or 'fwhm' other than the default 'vrad'.
    norm : string, optional
        The normalization; either "Scargle", "HorneBaliunas", or 
        "Cumming". Default is "HorneBaliunas".
    stats : boolean, optional
        Set True to obtain some statistical output (default is False).
  
  Attributes
    power : array
        The normalized power of the GLS.
    ofac : int
        The oversampling factor.
    hifac : float
        The maximum frequency.
  """

  def __init__(self, rv, ofac=6, hifac=1, period_range=None, time_of_offset=None, quantity='vrad',
               norm="HorneBaliunas", stats=False, full_output=False):

    if not gatspy_is_available:
      msg = red('ERROR: ') + 'gatspy is not available.'
      clogger.fatal(msg)
      return 

    self.name = 'Multiband Lomb-Scargle'
    try:
      self.star_name = get_star_name(rv)
    except AttributeError:
      self.star_name = ''

    self.power = None
    self.units = ''
    self.period_range = period_range
    self.ofac, self.hifac = ofac, hifac
    self.fullout = full_output
    self.t = rv.time
    self.th = rv.time - min(rv.time)

    if time_of_offset is None:
      # HARPS_upgrade
      self.time_of_offset = 57170
    else:
      assert isinstance(time_of_offset, (int, float)), '"time_of_offset" should be a number.'
      self.time_of_offset = time_of_offset

    msg = blue('INFO: ') + 'Time of offset is set to %d\n' % self.time_of_offset
    msg += blue('    : ') + 'Data before and after this time will be considered as different bands'
    clogger.info(msg)

    self.mask = self.t > self.time_of_offset

    if (~self.mask).all():
      msg = red('ERROR: ') + 'All observations are before %d. Multiband periodogram is not appropriate' % self.time_of_offset
      clogger.fatal(msg)
      raise ValueError

    if quantity == 'vrad':
      self.y = rv.vrad
      self.error = rv.error
      self.units = rv.units
    elif quantity == 'bis':
      self.y = rv.extras.bis_span
      self.error = 2. * rv.error #ones_like(self.y)
    elif quantity == 'fwhm':
      self.y = rv.extras.fwhm
      self.error = 2.35 * rv.error #ones_like(self.y)
    elif quantity == 'rhk':
      self.y = rv.extras.rhk
      self.error = rv.extras.sig_rhk
    elif quantity == 'contrast':
      self.y = rv.extras.contrast
      self.error = ones_like(self.y)
    elif quantity == 'resid':
      try:
        self.y = rv.fit['residuals']
      except TypeError:
        clogger.fatal('error!')
        return 
      self.error = rv.error

    self.norm = norm
    self._stats = stats

    self.__calcPeriodogram()
    
  def __calcPeriodogram(self):

    model = periodic.LombScargleMultiband(fit_period=True,
                                          optimizer_kwds={"quiet": True})

    if self.period_range is None:
      plow = max(0.5, 2*ediff1d(self.t).min()) # minimum meaningful period?
      phigh = self.t.ptp()
      model.optimizer.period_range = (plow, phigh)

    else:
      assert isinstance(self.period_range, tuple), 'period_range must be a tuple (plow, phigh)'
      assert len(self.period_range) == 2, 'period_range must be a tuple (plow, phigh)'
      plow = self.period_range[0]
      model.optimizer.period_range = self.period_range

    self.__buildFreq(plow=plow)
    periods = 1./self.freq

    model.fit(self.t, self.y, self.error, self.mask)
    self.power = model.periodogram(periods)

    offset = np.ediff1d(model.ymean_by_filt_)
    msg = yellow('RESULT: ') + 'Offset value: %f %s' % (offset[0], self.units)
    clogger.info(msg)

    self.model = model

    # self.N = len(self.y)
    # # An ad-hoc estimate of the number of independent frequencies 
    # # see discussion following Eq. (24)
    # self.M = (max(self.freq)-min(self.freq)) * (self.th.max() - self.th.min())

    # # Normalization:
    # if self.norm == "Scargle":
    #   popvar=raw_input('Input a priori known population variance:')
    #   self.power = self._upow/float(popvar)
    # if self.norm == "HorneBaliunas":
    #   self.power = (self.N-1.)/2.*self._upow
    # if self.norm == "Cumming":
    #   self.power = (self.N-3.)/2. * self._upow/(1.-max(self._upow))
    
    # Output statistics
    # if self._stats:
    #   self._output()
    # if self._showPlot:
    #   self._plot()

  def _normalize(self):
    if self.norm == "Scargle":
      popvar=raw_input('pyTiming::gls - Input a priori known population variance:')
      self.power = self._upow/float(popvar)
    if self.norm == "HorneBaliunas":
      self.power = (self.N-1.)/2.*self._upow
    if self.norm == "Cumming":
      self.power = (self.N-3.)/2. * self._upow/(1.-max(self._upow))

  def _normalize_value(self, a):
    # normalize the value (or array) a
    if self.norm == "Scargle":
      popvar=raw_input('pyTiming::gls - Input a priori known population variance:')
      return a/float(popvar)
    if self.norm == "HorneBaliunas":
      return (self.N-1.)/2.*a
    if self.norm == "Cumming":
      return (self.N-3.)/2. * a/(1.-max(np.atleast_1d(a)))    

  def __buildFreq(self, plow=0.5):
    """
      Build frequency array (`freq` attribute).
    """
    xdif = max(self.th)-min(self.th)
    # print xdif
    # nout = self.ofac * self.hifac * len(self.th)/2
    nout = int(xdif * self.ofac / plow)
    # print self.ofac, nout
    # sys.exit(0)
    self.freq = 1./(xdif) + arange(nout)/(self.ofac*xdif)




  
# This is the Bayesian Lomb-Scargle (BLS) periodogram as defined by
# Bretthorst (2001), Gregory (2005)
#####################################################################
class bls(PeriodogramBase):
  """
  Compute the Bayesian Lomb-Scargle (BLS) periodogram.

  This class implements the bayesian Lomb-Scargle periodogram as
  developed by Bretthorst (2000, 2001). This corresponds to the 
  bayesian expression for the posterior p(f|D,I).

  The constructor takes a RVSeries instance (i.e. a rv curve) as 
  first argument. The algorithm is slow due to the required frequency
  resolution; currently, only a Fortran implementation is available.
  There is an optional `freq` array, that can contain the 
  frequencies on which to calculate the periodogram. If not provided
  (....)
    
  Parameters
    rv : RVSeries
        The radial velocity curve or any object providing the attributes
        time, vrad and error which define the data.
    ofac : int
        Oversampling factor (default=6).
    hifac : float
        hifac * "average" Nyquist frequency is highest frequency for 
        which the periodogram will be calculated (default=40).
    freq : array, optional
        Contains the frequencies at which to calculate the periodogram.
        If not given, a frequency array will be automatically generated.
    quantity : string, optional
        For which quantity to calculate the periodogram. Possibilities are
        'bis', 'rhk', 'contrast' or 'fwhm' other than the default 'vrad'.    
    stats : boolean, optional
        Set True to obtain some statistical output (default is False).
  
  Attributes
    power : array
        The normalized power of the GLS.
    freq : array
        The frequency array.
  """

  def __init__(self, rv, ofac=6, hifac=40, freq=None, quantity='vrad', 
               stats=False, force_notrend=False):
    self.name = 'Bayesian Lomb-Scargle'

    self.power = None
    self.freq = freq
    self.ofac, self.hifac = ofac, hifac
    self.t = rv.time
    self.th = self.t - min(rv.time)
    if quantity == 'vrad':
      self.y = rv.vrad
      self.error = rv.error
    elif quantity == 'bis':
      self.y = rv.extras.bis_span
      self.error = 2. * rv.error #ones_like(self.y)
    elif quantity == 'fwhm':
      if not force_notrend and ask_yes_no('Should I remove a linear trend first? (y/N) ', False):
        m, b = polyfit(self.t, rv.extras.fwhm, 1)
        yp = polyval([m, b], self.t)
        self.y = rv.extras.fwhm - yp
      else:
        self.y = rv.extras.fwhm
      self.error = 2.35 * rv.error #ones_like(self.y)
    elif quantity == 'rhk':
      self.y = rv.extras.rhk
      self.error = rv.extras.sig_rhk
    elif quantity == 'contrast':
      self.y = rv.extras.contrast
      self.error = ones_like(self.y)

    # time span of observations
    self.Tspan = max(rv.time) - min(rv.time)
    # signal to noise ratio
    self.SNR = (0.5 * (max(self.y) - min(self.y)) ) / (sqrt(2.)*mean(self.error))

    
    self.__calcPeriodogram()
    
  def __calcPeriodogram(self):
    """ Compute the Bayesian Lomb-Scargle Periodogram. """

    # Build frequency array if not present
    if self.freq is None:
      plow = max(0.5, 2*ediff1d(self.t).min()) # minimum meaningful period?
      self.__buildFreq(plow=plow)
    # Circular frequencies
    omegas = 2.*pi * self.freq

    Pb, h2bar = blombscargle(self.t, self.y, self.freq)
    self.power = Pb 

  def __buildFreq(self, plow=0.5):
    """
      Build frequency array (`freq` attribute).
    """
    # the Bayesian version of the periodogram can produce arbitrarily
    # sharp spectral features. In general it's necessary to set the 
    # frequency resolution higher than for other periodograms. The 
    # guide is Eq. (13.4) of Gregory (2005)
    fhigh = 1./plow
    nyquist = 0.5*1./mean(self.t[1::]-self.t[0:-1])
    flow = 1./float(self.hifac) * nyquist
    # Eq. (13.4)
    df = 1. / (1.6 * self.SNR * self.Tspan * sqrt(len(self.y)))
    
    self.nf = int(0.1 * self.ofac * abs( (fhigh - flow) / df ))
    if self.nf > 1e6:
      self.nf = 1e6 # sensible limit to number of frequencies
    self.freq = linspace(flow, fhigh, self.nf) # could this be logspace?

  def prob(self, Pn):
    """
      Probability of obtaining the given power.
    
      Calculate the probability to obtain a power higher than
      `Pn` from the noise, which is assumed to be Gaussian.
      
      .. note:: Normalization
        (see [ZK09]_ for further details).

        - `Scargle`: 
        .. math::
          exp(-Pn)

        - `HorneBaliunas`: 
        .. math::
          \\left(1 - 2 \\times \\frac{Pn}{N-1} \\right)^{(N-3)/2}
        
        - `Cumming`: 
        .. math:: 
          \\left(1+2\\times \\frac{Pn}{N-3}\\right)^{-(N-3)/2}
      
      Parameters
        Pn : float
          Power threshold.
      
      Returns
        Probability : float
          The probability to obtain a power equal or
          higher than the threshold from the noise.
    """
    if self.norm=="Scargle": return exp(-Pn)
    if self.norm=="HorneBaliunas": return (1.-2.*Pn/(self.N-1.))**((self.N-3.)/2.)
    if self.norm=="Cumming": return (1.+2.*Pn/(self.N-3.))**(-(self.N-3.)/2.)

  def probInv(self, Prob):
    """
      Calculate minimum power for given probability.
    
      This function is the inverse of `Prob(Pn)`.
      Returns the minimum power for a given probability threshold Prob.
      
      Parameters
        Prob : float
          Probability threshold.
      
      Returns
        Power threshold : float
          The minimum power for the given
          false-alarm probability threshold.
    """
    if self.norm=="Scargle": return -log(Prob)
    if self.norm=="HorneBaliunas": return (self.N-1.)/2.*(1.-Prob**(2./(self.N-3.)))
    if self.norm=="Cumming": return (self.N-3.)/2.*(Prob**(-2./(self.N-3.))-1.) 


# This is the Bayesian Generalized Lomb-Scargle (BGLS) periodogram as 
# defined by Mortier et al (2014). 
#####################################################################
class bgls(PeriodogramBase):
  """
  Compute the Bayesian Generalized Lomb-Scargle (BGLS) periodogram.

  This class implements the Bayesian generalization of the 
  error-weighted Lomb-Scargle periodogram as developed by [M14]_.

  The constructor takes a RVSeries instance (i.e. a rv curve) as 
  first argument. As the algorithm is slow-ish, an implementation
  in Fortran is available if the keyword 'ext' is set to True in 
  the constructor or globally. 
  There is an optional `freq` array, that can contain the 
  frequencies on which to calculate the periodogram. If not provided
  (....)
    
  Parameters
    rv : RVSeries
        The radial velocity curve or any object providing the attributes
        time, vrad and error which define the data.
    ofac : int
        Oversampling factor (default=6).
    hifac : float
        hifac * "average" Nyquist frequency is highest frequency for 
        which the periodogram will be calculated (default=1).
    freq : array, optional
        Contains the frequencies at which to calculate the periodogram.
        If not given, a frequency array will be automatically generated.
    quantity : string, optional
        For which quantity to calculate the periodogram. Possibilities are
        'bis', 'rhk', 'contrast' or 'fwhm' other than the default 'vrad'.
    norm : string, optional
        The normalization; either "Scargle", "HorneBaliunas", or 
        "Cumming". Default is "HorneBaliunas".
    stats : boolean, optional
        Set True to obtain some statistical output (default is False).
    ext : boolean, optional
        Use Fortran extension in the calculation (default is True)
  
  Attributes
    power : array
        The normalized power of the GLS.
    freq : array
        The frequency array.
    ofac : int
        The oversampling factor.
    hifac : float
        The maximum frequency.
    norm : string
        The normalization used.
  """

  def __init__(self, rv, ofac=6, hifac=1, freq=None, quantity='vrad', stats=False, ext=False):
    self.name = 'Bayesian Generalized Lomb-Scargle'

    self.power = None
    self.freq = freq
    self.ofac, self.hifac = ofac, hifac
    self.t = rv.time
    self.th = rv.time - min(rv.time)
    if quantity == 'vrad':
      self.y = rv.vrad
      self.error = rv.error
    elif quantity == 'bis':
      self.y = rv.extras.bis_span
      self.error = 2. * rv.error #ones_like(self.y)
    elif quantity == 'fwhm':
      if ask_yes_no('Should I remove a linear trend first? (y/N) ', False):
        m, b = polyfit(self.t, rv.extras.fwhm, 1)
        yp = polyval([m, b], self.t)
        self.y = rv.extras.fwhm - yp
      else:
        self.y = rv.extras.fwhm
      self.error = 2.35 * rv.error #ones_like(self.y)
    elif quantity == 'rhk':
      self.y = rv.extras.rhk
      self.error = rv.extras.sig_rhk
    elif quantity == 'contrast':
      self.y = rv.extras.contrast
      self.error = ones_like(self.y)
    elif quantity == 'resid':
      try:
        self.y = rv.fit['residuals']
      except TypeError:
        clogger.fatal('error!')
        return 
      self.error = rv.error

    # Check and assign normalization
    self.label = {'title': 'Generalized Lomb Periodogram',\
                  'xlabel': 'Frequency',\
                  'ylabel': "Posterior pdf"}
    
    self._stats = stats
    if ext: 
      self.__calcPeriodogramFast()
    else: 
      self.__calcPeriodogram()
    
  def __calcPeriodogram(self):
    """ Compute the BGLS. Notation very close to [M14]_ """
    import mpmath

    # Build frequency array if not present
    if self.freq is None:
      plow = max(0.5, 2*ediff1d(self.t).min()) # minimum meaningful period?
      self.__buildFreq(plow=plow)
    # Circular frequencies
    omegas = 2.*pi * self.freq

    err2 = self.error*self.error
    w = 1./err2
    W = sum(w)


    bigY = sum(w*self.y)  # Eq. (10)

    p = []
    constants = []
    exponents = []

    for i, omega in enumerate(omegas):
      theta = 0.5 * arctan2(sum(w*sin(2.*omega*self.t)), sum(w*cos(2.*omega*self.t)))
      x = omega*self.t - theta
      cosx = cos(x)
      sinx = sin(x)
      wcosx = w*cosx
      wsinx = w*sinx

      C = sum(wcosx)
      S = sum(wsinx)

      YCh = sum(self.y*wcosx)
      YSh = sum(self.y*wsinx)
      CCh = sum(wcosx*cosx)
      SSh = sum(wsinx*sinx)

      if (CCh != 0 and SSh != 0):
        K = (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2.*CCh*SSh)
        L = (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
        M = (YCh*YCh*SSh + YSh*YSh*CCh)/(2.*CCh*SSh)
        constants.append(1./sqrt(CCh*SSh*abs(K)))
      elif (CCh == 0):
        K = (S*S - W*SSh)/(2.*SSh)
        L = (bigY*SSh - S*YSh)/(SSh)
        M = (YSh*YSh)/(2.*SSh)
        constants.append(1./sqrt(SSh*abs(K)))
      elif (SSh == 0):
        K = (C*C - W*CCh)/(2.*CCh)
        L = (bigY*CCh - C*YCh)/(CCh)
        M = (YCh*YCh)/(2.*CCh)
        constants.append(1./sqrt(CCh*abs(K)))
      if K > 0:
        raise RuntimeError('K is positive. This should not happen.')

      exponents.append(M - L*L/(4.*K))

    constants = array(constants)
    exponents = array(exponents)

    logp = log10(constants) + (exponents * log10(exp(1.)))

    p = [10**mpmath.mpf(x) for x in logp]

    msg = yellow('RESULT: ') + 'Maximum of posterior distribution - %s' % mpmath.nstr(max(p))
    clogger.info(msg)
    msg = blue('INFO: ') + 'Normalizing by this value'
    clogger.info(msg)


    p = array(p) / max(p)  # normalize

    p[p < (float_info.min * 10)] = 0
    self._upow = array([float(pp) for pp in p])
    self.power = self._upow

    self.N = len(self.y)
    # An ad-hoc estimate of the number of independent frequencies 
    # see discussion following [ZK09]_ Eq. (24)
    self.M = (max(self.freq)-min(self.freq)) * (self.th.max() - self.th.min())

    # Output statistics
    if self._stats:
      self._output()
    # if self._showPlot:
    #   self._plot()

  def __calcPeriodogramFast(self):
    raise NotImplementedError('There is no fast version of BGLS yet.')


  def __buildFreq(self, plow=0.5):
    """
      Build frequency array (`freq` attribute).
    """
    phigh = self.th[-1]
    n_steps = int(self.ofac * self.th.size * (1./plow - 1./phigh))
    # print n_steps
    # xdif = max(self.th)-min(self.th)
    ##### nout = self.ofac * self.hifac * len(self.th)/2
    # nout = int(xdif * self.ofac / plow)
    # self.freq = 1./(xdif) + arange(nout)/(self.ofac*xdif)
    self.freq = np.linspace(1./phigh, 1./plow, n_steps)

  def prob(self, Pn):
    raise NotImplementedError('Not implemented for BGLS yet.')

  def probInv(self, Prob):
    raise NotImplementedError('Not implemented for BGLS yet.')


# This is the Hoeffding-test periodicity metric based upon 
# Hoeffding 1948, Ann. Math. Stat. 19, 293
#####################################################################
class hoeffding(PeriodogramBase):
  """
  Calculate The Hoeffding-test "periodogram".

  This class implements the Hoeffding-test periodicity metric based
  upon Hoeffding 1948, Ann. Math. Stat. 19, 293 and as developed by
  Shay Zucker 2015, arXiv:1503.01734

  The constructor takes a RVSeries instance (i.e. a rv curve) as 
  first argument.
  There is an optional `freq` array, that can contain the 
  frequencies on which to calculate the periodogram. If not provided
  (....)
    
  Parameters
    rv : RVSeries
        The radial velocity curve or any object providing the attributes
        time, vrad and error which define the data.
    ofac : int
        Oversampling factor (default=6).
    hifac : float
        hifac * "average" Nyquist frequency is highest frequency for 
        which the periodogram will be calculated (default=1).
    freq : array, optional
        Contains the frequencies at which to calculate the periodogram.
        If not given, a frequency array will be automatically generated.
    quantity : string, optional
        For which quantity to calculate the periodogram. Possibilities are
        'bis', 'rhk', 'contrast' or 'fwhm' other than the default 'vrad'.
    ext : boolean, optional
        Use Fortran extension in the calculation (default is False)
  
  Attributes
    power : array
        The normalized power of the GLS.
    freq : array
        The frequency array.
    ofac : int
        The oversampling factor.
    hifac : float
        The maximum frequency.
    norm : string
        The normalization used.
  """

  def __init__(self, rv, ofac=4, hifac=1, freq=None, quantity='vrad', ext=False):
    self.name = 'Hoeffding'

    self.power = None
    self.freq = freq
    self.ofac, self.hifac = ofac, hifac
    self.t = rv.time
    self.th = rv.time - min(rv.time)
    if quantity == 'vrad':
      self.y = rv.vrad
      self.error = rv.error
    elif quantity == 'bis':
      self.y = rv.extras.bis_span
      self.error = 2. * rv.error #ones_like(self.y)
    elif quantity == 'fwhm':
      # if not force_notrend and ask_yes_no('Should I remove a linear trend first? (y/N) ', False):
      #   m, b = polyfit(self.t, rv.extras.fwhm, 1)
      #   yp = polyval([m, b], self.t)
      #   self.y = rv.extras.fwhm - yp
      # else:
      self.y = rv.extras.fwhm
      self.error = 2.35 * rv.error #ones_like(self.y)
    elif quantity == 'rhk':
      self.y = rv.extras.rhk
      self.error = rv.extras.sig_rhk
    elif quantity == 'contrast':
      self.y = rv.extras.contrast
      self.error = ones_like(self.y)
    elif quantity == 'resid':
      try:
        self.y = rv.fit['residuals']
      except TypeError:
        clogger.fatal('error!')
        return 
      self.error = rv.error

    self.label = {'title': self.name + ' Periodogram'}
    
    if ext: 
      self.__calcPeriodogramFast()
    else: 
      self.__calcPeriodogram()
    
  def __calcPeriodogram(self):
    """ Compute the Hoeffding-test """
    # Build frequency array if not present
    if self.freq is None:
      plow = 0.5#max(0.5, 2*ediff1d(self.t).min()) # minimum meaningful period?
      self.__buildFreq(plow=plow)

    N = self.t.size
    Nfreqs = self.freq.size

    # Temporaries
    t = self.t
    y = self.y


    # Allocate the output vector (to save memory handling time)
    D = zeros_like(self.freq)

    # Allocate area for auxiliary vectors in advance
    a = zeros_like(t) #zeros(N,1);
    b = zeros_like(t) #zeros(N,1);
    c = zeros_like(t) #zeros(N,1);

    periods = 1./self.freq

    # Convert the original values to ranks, using Scipy's rankdata
    rankx = rankdata(y)  #tiedrank(x);

    # Serially enumerate the periods (frequencies)
    for ifreq in range(Nfreqs):
      # Phase folding
      phases = np.mod(t, periods[ifreq])

      iphases = np.argsort(phases)
      sortedphases = phases[iphases]

      # A special treatment for the wrap-around.
      iphases = np.append(iphases, iphases[0])

      # Calculate Hoeffding test statistic (D)
      # --------------------------------------
      # a = 'current' sample rank
      # b = 'next' sample rank (in terms of the phase ordering)
      a = rankx[iphases[:-1]] - 1 # rankx(iphases(1:(end-1)))-1;
      b = rankx[iphases[1:]] - 1 # rankx(iphases(2:end))-1;
      
      # Calculate the bivariate rank
      for ii in range(N):
        c[ii] = sum((y[iphases[:-1]] < y[iphases[ii]]) & (y[iphases[1:]] < y[iphases[ii+1]]))

      # Calculate the A, B and C statistics in Hoeffding's formulae
      A = sum(a*(a-1)*b*(b-1)) # sum(a.*(a-1).*b.*(b-1));
      B = sum((a-1)*(b-1)*c) # sum((a-1).*(b-1).*c);
      C = sum(c*(c-1)) # sum(c.*(c-1));
      
      D[ifreq] = A-2*(N-2)*B+(N-2)*(N-3)*C

      # Save the current phase ordering, assuming next period's ordering 
      # will not be too much different, thus saving some time.
      t = t[iphases[:-1]] # t(iphases(1:end-1));
      y = y[iphases[:-1]] # x(iphases(1:end-1));
      rankx = rankx[iphases[:-1]] # rankx(iphases(1:end-1));
    
    # Apply Hoeffding's normalization
    D = D/(N*(N-1)*(N-2)*(N-3)*(N-4))

    self._a = a
    self._b = b
    self._c = c

    self.power = D


  def __calcPeriodogramFast(self):
    """ """
    pass

  def _normalize(self):
    pass

  def _normalize_value(self, a):
    pass

  def __buildFreq(self, plow=0.5):
    """ Build frequency array (`freq` attribute) """
    xdif = max(self.th)-min(self.th)
    # print xdif
    # nout = self.ofac * self.hifac * len(self.th)/2
    nout = int(xdif * self.ofac / plow)
    # print nout
    # sys.exit(0)
    self.freq = 1./(xdif) + arange(nout)/(self.ofac*xdif)

  def prob(self, Pn):
    pass

  def probInv(self, Prob):
    pass


# This is the ...... periodogram etc etc etc
#
#
#####################################################################
class SpectralWindow(PeriodogramBase):
  """ Calculate the spectral window function and its phase angles.
  See Eq. (1) of Dawson & Fabrycky (2010).

    Parameters
      freq : float, array
        Frequencies at which to calculate the spectral window.
      time : float, array
        Times of the data points.
    
    Returns
      amp : float, array, len(freq)
        Amplitude of the spectral window function.
      phase : float, array, len(freq)
        Phase angles of the spectral window function.
  """

  def __init__(self, freq, time):
    self.name = 'SpectralWindow'
    self.freq = freq
    self.time = time
    self.amp = None
    self.phase = None
    self._calcWindowFunction()
    self._plot(dials=False)


  def _calcWindowFunction(self, new=False):
    n = self.time.size

    if not new:
      self.W = np.fromiter((np.sum(np.exp(-2.j*pi*f*self.time))/n for f in self.freq), np.complex_, self.freq.size)
      # self.W = array([sum([exp(-2.j*pi*f*t) for t in self.time])/float(n) for f in self.freq])
      self.amp = np.absolute(self.W)
      self.power = self.amp
      self.phase = np.arctan2(self.W.imag, self.W.real)
    else:
      self.W2 = np.fromiter((np.sum(np.exp(-2.j*pi*f*self.time))/n for f in self.freq2), np.complex_, self.freq2.size)
      self.amp2 = np.absolute(self.W2)
      self.power2 = self.amp2
      self.phase2 = np.arctan2(self.W2.imag, self.W2.real)


  def _plot(self, dials=False, ndials=3):
    """
      Plot the spectral window function as a function of period.
    """
    self.fig = plt.figure()
    self.ax = self.fig.add_subplot(1,1,1)
    # self.ax.set_title("Spectral Window Function")
    self.ax.set_xlabel("Period")
    self.ax.set_ylabel("Window function")
    self.ax.semilogx(1./self.freq, self.amp, 'b-')

    # Trying to get labels on top representing frequency - does not work!
    # self.ax2 = self.ax.twiny()
    # def tick_function(X):
    #   V = 1./X
    #   return ["%.e" % z for z in V]

    # bottom_tick_loc = self.ax.get_xticks()
    # self.ax2.set_xticks(bottom_tick_loc)
    # self.ax2.set_xticklabels(tick_function(bottom_tick_loc))
    # self.ax2.set_xlabel(r"Frequency [day$^{-1}$]")
    # # self.ax2.get_xaxis().get_major_formatter().set_scientific(True)

    if dials:
      fmax1, fmax2, fmax3 = self.get_peaks(n=3)
      max_amp = max(self.amp)

      self.ax.semilogx(1./fmax1,max_amp+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
      self.ax.semilogx(1./fmax3,max_amp+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)
      self.ax.semilogx(1./fmax2,max_amp+0.1,marker = '$\circ$',markersize=10,c='k',mew=0.3)

      ph1 = ph2 = ph3 = 0.3

      self.ax.semilogx([1./fmax1,1./fmax1+0.025*cos(ph1)],[max_amp+0.1,max_amp+0.1+0.025*sin(ph1)],'k-',lw=1)
      self.ax.semilogx([1./fmax2,1./fmax2+0.025*cos(ph2)],[max_amp+0.1,max_amp+0.1+0.025*sin(ph2)],'k-',lw=1)
      self.ax.semilogx([1./fmax3,1./fmax3+0.025*cos(ph3)],[max_amp+0.1,max_amp+0.1+0.025*sin(ph3)],'k-',lw=1)


    # self.fig.tight_layout()
    plt.show()


  def _plot_freq(self, dials=False, ndials=3):
    """
      Plot the spectral window function as a function of frequency.
    """
    # we need different trial frequencies for this plot
    self.freq2 = np.linspace(1e-5, 4.2, 1e4)
    self._calcWindowFunction(new=True)
    self.fig = plt.figure()
    self.ax = self.fig.add_subplot(1,1,1)
    # self.ax.set_title("Spectral Window Function")
    self.ax.set_xlabel("Frequency")
    self.ax.set_ylabel("Window function")
    self.ax.plot(self.freq2, self.amp2, 'b-')

    # Trying to get labels on top representing frequency - does not work!
    # self.ax2 = self.ax.twiny()
    # def tick_function(X):
    #   V = 1./X
    #   return ["%.e" % z for z in V]

    # bottom_tick_loc = self.ax.get_xticks()
    # self.ax2.set_xticks(bottom_tick_loc)
    # self.ax2.set_xticklabels(tick_function(bottom_tick_loc))
    # self.ax2.set_xlabel(r"Frequency [day$^{-1}$]")
    # # self.ax2.get_xaxis().get_major_formatter().set_scientific(True)

    if dials:
      fmax1, fmax2, fmax3 = self.get_peaks(n=3)
      max_amp = max(self.amp2)

      self.ax.semilogx(fmax1, max_amp+0.1, marker='$\circ$', markersize=10, c='k', mew=0.3)
      self.ax.semilogx(fmax3, max_amp+0.1, marker='$\circ$', markersize=10, c='k', mew=0.3)
      self.ax.semilogx(fmax2, max_amp+0.1, marker='$\circ$', markersize=10, c='k', mew=0.3)

      ph1 = ph2 = ph3 = 0.3

      self.ax.semilogx([fmax1, fmax1+0.025*cos(ph1)], [max_amp+0.1, max_amp+0.1+0.025*sin(ph1)], 'k-', lw=1)
      self.ax.semilogx([fmax2, fmax2+0.025*cos(ph2)], [max_amp+0.1, max_amp+0.1+0.025*sin(ph2)], 'k-', lw=1)
      self.ax.semilogx([fmax3, fmax3+0.025*cos(ph3)], [max_amp+0.1, max_amp+0.1+0.025*sin(ph3)], 'k-', lw=1)


    # self.fig.tight_layout()
    plt.show()


#####################################################################
### References 
#####################################################################
# [ZK09] - Zechmeister, M. & Kürster, M., A&A 496, 577-584, 2009