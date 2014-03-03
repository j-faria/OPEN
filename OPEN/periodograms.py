# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
from classes import PeriodogramBase
from numpy import *
from numpy.fft import *  
import cmath
import matplotlib.pyplot as plt 

from .logger import clogger, logging

from ext.blombscargle import blombscargle
from ext.glombscargle import glombscargle



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
  def __init__(self, rv, ofac=6, hifac=6):
    self.name = 'Lomb-Scargle - Press & Rybicki algorithm'

    self.power = None
    self.freq = None
    self.ofac = ofac
    self.hifac = hifac
    self.t = rv.time
    self.y = rv.vrad
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
  ----------
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
        'bis' or 'fwhm' other than the default 'vrad'.
    norm : string, optional
        The normalization; either "Scargle", "HorneBaliunas", or 
        "Cumming". Default is "HorneBaliunas".
    stats : boolean, optional
        Set True to obtain some statistical output (default is False).
    ext : boolean, optional
        Use Fortran extension in the calculation (default is True)
  
  Attributes
  ----------
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
               norm="HorneBaliunas", stats=False, ext=True):
    self.name = 'Generalized Lomb-Scargle'

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
      self.error = ones_like(self.y)
    elif quantity == 'fwhm':
      self.y = rv.extras.fwhm
      self.error = 2.35 * rv.error #ones_like(self.y)
    elif quantity == 'rhk':
      self.y = rv.extras.rhk
      self.error = rv.extras.sig_rhk
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
      self.__buildFreq()
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
      popvar=raw_input('pyTiming::gls - Input a priori known population variance:')
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

    # Build frequency array if not present
    if self.freq is None:
      self.__buildFreq()
    # Circular frequencies
    omegas = 2.*pi * self.freq

    # unnormalized power and an estimate of the number of independent frequencies 
    self._upow, self.M = glombscargle(self.t, self.y, self.error, omegas)

    self.N = len(self.y)
    # Normalization:
    if self.norm == "Scargle":
      popvar=raw_input('pyTiming::gls - Input a priori known population variance:')
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

  def __buildFreq(self, plow=0.5):
    """
      Build frequency array (`freq` attribute).
    """
    xdif = max(self.th)-min(self.th)
    # print xdif
    # nout = self.ofac * self.hifac * len(self.th)/2
    nout = int(xdif * self.ofac / plow)
    # print nout
    # sys.exit(0)
    self.freq = 1./(xdif) + arange(nout)/(self.ofac*xdif)

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
      ----------
      Pn : float
          Power threshold.
      
      Returns
      -------
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
      ----------
      Prob : float
          Probability threshold.
      
      Returns
      -------
      Power threshold : float
          The minimum power for the given
          false-alarm probability threshold.
    """
    if self.norm=="Scargle": return -log(Prob)
    if self.norm=="HorneBaliunas": return (self.N-1.)/2.*(1.-Prob**(2./(self.N-3.)))
    if self.norm=="Cumming": return (self.N-3.)/2.*(Prob**(-2./(self.N-3.))-1.)

  
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
  ----------
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
        'bis' or 'fwhm' other than the default 'vrad'.    
    stats : boolean, optional
        Set True to obtain some statistical output (default is False).
  
  Attributes
  ----------
    power : array
        The normalized power of the GLS.
    freq : array
        The frequency array.
  """

  def __init__(self, rv, ofac=6, hifac=40, freq=None, quantity='vrad', stats=False):
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
      self.error = mean(rv.error) * ones_like(self.y)
    elif quantity == 'fwhm':
      self.y = rv.extras.fwhm
      self.error = mean(rv.error) * ones_like(self.y)

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
    
    self.nf = int(0.1 * abs( (fhigh - flow) / df ))
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
      ----------
      Pn : float
          Power threshold.
      
      Returns
      -------
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
      ----------
      Prob : float
          Probability threshold.
      
      Returns
      -------
      Power threshold : float
          The minimum power for the given
          false-alarm probability threshold.
    """
    if self.norm=="Scargle": return -log(Prob)
    if self.norm=="HorneBaliunas": return (self.N-1.)/2.*(1.-Prob**(2./(self.N-3.)))
    if self.norm=="Cumming": return (self.N-3.)/2.*(Prob**(-2./(self.N-3.))-1.) 


# This is the ...... periodogram etc etc etc
#
#
#####################################################################
class SpectralWindow(PeriodogramBase):
  """ Calculate the spectral window function and its phase angles.
  See Eq. (1) of Dawson & Fabrycky (2010).

    Parameters
    ----------
    freq : float, array
        Frequencies at which to calculate the spectral window.
    time : float, array
        Times of the data points.
    
    Returns
    -------
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
    self._plot(dials=True)


  def _calcWindowFunction(self):
    n = len(self.time)
    W = [sum([cmath.exp(-2.j*pi*f*t) for t in self.time])/float(n) for f in self.freq]

    self.amp = [sqrt(w.real*w.real + w.imag*w.imag) for w in W]
    self.phase = [arctan2(w.imag, w.real) for w in W]


  def _plot(self, dials=False, ndials=3):
    """
      Plot the spectral window function as a function of frequency.
    """
    self.fig = plt.figure()
    self.ax = self.fig.add_subplot(1,1,1)
    # self.ax.set_title("Spectral Window Function")
    self.ax.set_xlabel("Period")
    self.ax.set_ylabel("Power")
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


    self.fig.tight_layout()
    plt.show()
    # p = pg.plot(1./self.freq, self.power, title="Periodogram")
    # p.plotItem.setLogMode(1./self.freq, self.power)
    # pg.QtGui.QApplication.exec_()



#####################################################################
### References 
#####################################################################
# [ZK09] - Zechmeister, M. & Kürster, M., A&A 496, 577-584, 2009