/* 
Adapted from code by Joel Hartman by João Faria
http://www.astro.princeton.edu/~jhartman/vartools.html 

Further adapted by Annelies Mortier
*/
#include "Python.h"
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>

#define DFTTINY 1.0e-15
#define BEAMTOL 1.0e-5
#define MIN_(A,B) ((A) < (B) ? (A) : (B))
#define MAX_(A,B) ((A) > (B) ? (A) : (B))

void nrerror(char error_text[])
{
/*
Numerical Recipes standard error handler
*/
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


void FDFT(double *x, double *y, double *sig, int N, double df, double Nf, double *f_r, double *f_i)
/* This computes the fast DFT of a signal y with uneven time sampling x and data errors sig
using the algorithm from Kurtz D.W. 1985, MNRAS, 213, 773. 
The number of points in x is N, the number of output frequencies is 2*Nf + 1 with a step-size of df. 
The real part of the DFT is stored in f_r, the imaginary part in f_i, 

output = f_r + if_i = 1/(NW) sum[wj*(yj - meany)*exp(-2pi i nu (xj - meanx) )]

where wj = 1/sigj^2 and W is the sum of wj

Remember: exp(i theta) = cos(theta) + i sin(theta)
*/
{
  int i, j, Ntot;
  double avx, avy, c, s, yval, wj, W;
  double xf0, xdf, sdf, cdf, s0, c0, xt, f0;

  /* First get the average x and y values */
  for(i=0, avx = 0., avy = 0.;i<N;i++)
    {
      avx += x[i];
      avy += y[i];
    }
  avx = avx/N;
  avy = avy/N;

  /* Define the weights sum */
  for(i=0,W=0;i<N;i++)
    {
      W += 1.0/(sig[i]*sig[i]);
    }

  /* Define the amount of freqs and initialize result arrays */ 
  Ntot = 2*Nf + 1;
  for(i=0;i<Ntot;i++)
    {
      f_r[i] = 0.;
      f_i[i] = 0.;
    }

  /* for each point in the data, get the corresponding output value */
  for(j=0;j<N;j++)
    {
      f0 = -Nf*df; // startfreq
      xt = 2.0*M_PI * (x[j] - avx);
      xf0 = xt*f0;
      xdf = xt*df;
      s0 = sin(xf0);
      c0 = cos(xf0);
      sdf = sin(xdf);
      cdf = cos(xdf);
      yval = y[j] - avy;
      wj = 1.0/(sig[j]*sig[j]);
      for(i=0;i<Ntot;i++)
	{
	  c = c0*cdf - s0*sdf;
	  s = s0*cdf + c0*sdf;
	  f_r[i] += wj*yval*c;
	  f_i[i] += wj*yval*s;
	  c0 = c;
	  s0 = s;
	}
    }
  for(i=0;i<Ntot;i++)
    {
      f_r[i] /= (double) N;
      f_i[i] /= (double) N;
      f_r[i] /= (double) W;
      f_i[i] /= (double) W;
    }
}

void FDFT_Wfunc(double *x, int N, double df, double Nf, double *W_r, double *W_i)
/* Same as above, this time for computing the window spectrum (no y) 

Wnu = W_r + i W_i = 1/N  sum[exp(-2pi i nu (xj - meanx) )]

Remember: exp(i theta) = cos(theta) + i sin(theta)
*/
{
  int i, j, Ntot;
  double avx, c, s;
  double xf0, xdf, sdf, cdf, s0, c0, xt, f0;

  /* First get the average x values */
  for(i=0, avx = 0.;i<N;i++)
    {
      avx += x[i];
    }
  avx = avx/N;

  Ntot = 2*Nf + 1;
  for(i=0;i<Ntot;i++)
    {
      W_r[i] = 0.;
      W_i[i] = 0.;
    }
  
  for(j=0;j<N;j++)
    {
      f0 = -Nf*df;
      xt = 2.0*M_PI * (x[j] - avx);
      xf0 = xt*f0;
      xdf = xt*df;
      s0 = sin(xf0);
      c0 = cos(xf0);
      sdf = sin(xdf);
      cdf = cos(xdf);
      for(i=0;i<Ntot;i++)
	{
	  c = c0*cdf - s0*sdf;
	  s = s0*cdf + c0*sdf;
	  W_r[i] += c;
	  W_i[i] += s;
	  c0 = c;
	  s0 = s;
	}
    }
  for(i=0;i<Ntot;i++)
    {
      W_r[i] /= (double) N;
      W_i[i] /= (double) N;
    }
}

void getamplitude(double *f_r, double *f_i, double *W_r, double *W_i, int Nf, int peakfindx, double *a_r, double *a_i)
/* Computes equation 25 of Roberts et al. 1987, alpha!! of a spectrum f and a windowfunction W

alpha(f,W,nu) = (fnu - fnuconjugate * W(2nu)) / (1 - amp(W(2nu))^2)

The real part of the amplitude will be stored in a_r, the imaginary part in a_i, 
Nf is the number of frequency points in the f spectrum.
peakfindx is the index pointing to the frequency we want to get the amplitude of
*/
{
  double denom;
  int i_f, i_w;
  
  i_f = peakfindx;
  i_w = 2*peakfindx;
 
  denom = 1. - W_r[i_w]*W_r[i_w] - W_i[i_w]*W_i[i_w];
  if(denom < DFTTINY)
    {
      *a_r = 0.;
      *a_i = 0.;
    }
  else
    {
      *a_r = (f_r[i_f] - f_r[i_f]*W_r[i_w] - f_i[i_f]*W_i[i_w])/denom;
      *a_i = (f_i[i_f] + f_i[i_f]*W_r[i_w] - f_r[i_f]*W_i[i_w])/denom;
    }
}

void doclean(double *R_r, double *R_i, double *W_r, double *W_i, int Nf, double *C_r, double *C_i, double gainval, double SNlimit, double df, int *iter)
/* This runs the clean algorithm on the dirty spectrum R with window spectrum W. 
The output cleaned spectrum is in C, with R being replaced by the final residual spectrum. 
Note that we do not find the clean beam or add the residual spectrum to the clean spectrum 
at this stage - that will be done later. 
The iteration stops if the peak is less than SNlimit times the noise OR if the residual
spectrum does not change enough anymore (defined by DFTTINY)
*/
{
  int i, j, Nftot;
  int peaki;
  double peakval, rmsval, ave1, ave2, val, a_r, a_i, lastsum, thissum;

  Nftot = 2*Nf + 1;
  /* Initialize the clean spectrum to 0. */
  for(i=0;i<Nftot;i++)
    {
      C_r[i] = 0.;
      C_i[i] = 0.;
    }

  lastsum = 0.;
  (*iter) = 0;
  while(1)
    {

      (*iter)++;

      /* First find the peak in the residual power spectrum (searching only positive frequencies), also compute the RMS of the power spectrum */
      for(i=Nf+1, peaki=-1, peakval = -1., ave1 = 0., ave2 = 0.;i<Nftot;i++)
	{
	  val = R_r[i]*R_r[i] + R_i[i]*R_i[i];
	  if(val > peakval)
	    {
	      peakval = val;
	      peaki = i;
	    }
	  ave1 += val;
	  ave2 += val*val;
	}
      ave1 = ave1 / (double) (Nf - 1);
      ave2 = ave2 / (double) (Nf - 1);
      rmsval = sqrt(ave2 - (ave1*ave1));

      /* Stop the iteration if the peak is less than SNlimit times the noise */
      if(peakval - ave1 < SNlimit * rmsval) break;

      /* Get the amplitude of the peak */
      getamplitude(R_r, R_i, W_r, W_i, Nf, peaki, &a_r, &a_i);
      
      peaki -= Nf;
      /* Calculate the new residual spectrum */
      for(j=0;j<Nftot;j++)
	{
	  R_r[j] -= gainval*(a_r*W_r[Nf+j-peaki] - a_i*W_i[Nf+j-peaki] + a_r*W_r[Nf+j+peaki] + a_i*W_i[Nf+j+peaki]);
	  R_i[j] -= gainval*(a_i*W_r[Nf+j-peaki] + a_r*W_i[Nf+j-peaki] - a_i*W_r[Nf+j+peaki] + a_r*W_i[Nf+j+peaki]);
	}

      thissum = lastsum;

      /* Update the clean spectrum */
      thissum -= (C_r[Nf + peaki]*C_r[Nf + peaki] + C_i[Nf + peaki]*C_i[Nf + peaki]);
      C_r[Nf + peaki] += gainval*a_r;
      C_i[Nf + peaki] += gainval*a_i;
      C_r[Nf - peaki] += gainval*a_r;
      C_i[Nf - peaki] -= gainval*a_i;
      thissum += (C_r[Nf + peaki]*C_r[Nf + peaki] + C_i[Nf + peaki]*C_i[Nf + peaki]);

      /* Break if the clean spectrum isn't changing appreciably any more */
      if(thissum < DFTTINY && lastsum < DFTTINY)
	break;
      if(thissum - lastsum < DFTTINY)
	break;
      if(lastsum > 0 ? (thissum < DFTTINY*lastsum) : 0)
	break;
      lastsum=thissum;
    }
}

void GetCleanBeam(double *W_r, double *W_i, int Nf, double **B_r, int *Nb)
/* This function fits a gaussian to the central peak of the window function, 
it will allocate memory for the beam and automatically determine an appropriate size for it
(step 4 Roberts et al.)
*/
{
  int i, j, Nfit, Nftot, x, nincrease;
  double centval, lastval, val, sigval, x2, f2, f3, val1, val2;

  /* First determine an appropriate size for fitting the beam */
  centval = sqrt(W_r[2*Nf]*W_r[2*Nf] + W_i[2*Nf]*W_i[2*Nf]);
  Nftot = 4*Nf + 1;
  lastval = centval;
  nincrease = 0;
  for(i=2*Nf+1,j=0;i<Nftot;i++,j++)
    {
      val = sqrt(W_r[i]*W_r[i] + W_i[i]*W_i[i]);
      if(val < BEAMTOL * centval || (val > lastval && nincrease > 2))
	break;
      if(val > lastval)
	nincrease++;
      else
	nincrease = 0;
      lastval = val;
    }
  Nfit = j;
  
  val1 = 0.;
  val2 = 0.;
  for(i=0,j=2*Nf-(Nfit),x=-(Nfit);i<(2*(Nfit) + 1);i++,j++,x++)
    {
      x2 = (double) x*x;
      val = sqrt(W_r[j]*W_r[j] + W_i[j]*W_i[j]);
      f2 = val*val;
      f3 = log(val)*f2;
      val1 -= x2*f3;
      val2 += x2*x2*f2;
    }
  sigval = val1 / val2;

  (*Nb) = ceil(sqrt(11.5/sigval));

  if(((*B_r) = (double *) malloc((2*(*Nb) + 1)* sizeof(double))) == NULL)
    nrerror("Memory allocation problem");

  for(i=0,x=-(*Nb);i<(2*(*Nb) + 1);i++,x++)
    {
      x2 = (double) -x*x*sigval;
      (*B_r)[i] = exp(x2);
    }
}

void Convolve_CleanBeam(double *C_r, double *C_i, double *R_r, double *R_i, double *S_r, double *S_i, int Nf, int Nb, double *B_r)
/*
Convolve the clean beam with the cleaned components and add the final residual spectrum
(step 5 Roberts et al.)
*/
{
  int j, k,l,Nftot;
  Nftot = 2*Nf + 1;
  for(j=0;j<Nftot;j++)
    {
      S_r[j] = R_r[j];
      S_i[j] = R_i[j];
      for(k=MAX_((j-Nb),0);k<MIN_((j+Nb),Nftot);k++, l++)
	{
	  l = Nb + j - k;
	  S_r[j] += C_r[k]*B_r[l];
	  S_i[j] += C_i[k]*B_r[l];
	}
    }
}

static PyObject *clean(PyObject *self, PyObject *args)
{
/* The python connected code to run the CLEAN algorithm */

    // input
    PyObject *t_obj;	// time
    PyObject *vrad_obj;	// RV
    PyObject *sig_obj;	// errors

    double plow;	// Lowest period to sample

    /* Argument handling */
    if ( !PyArg_ParseTuple(args, "O!O!O!d",
            &PyArray_Type, &t_obj,
            &PyArray_Type, &vrad_obj,
            &PyArray_Type, &sig_obj,
            &plow)
       )
    {
        return NULL;
    }

    /* Interpret the input objects as numpy arrays. */
    PyObject *t_array = PyArray_FROM_OTF(t_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *vrad_array = PyArray_FROM_OTF(vrad_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *sig_array = PyArray_FROM_OTF(sig_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    /* If that didn't work, throw an exception. */
    if (t_array == NULL || vrad_array == NULL || sig_array == NULL) {
      Py_XDECREF(t_array);
      Py_XDECREF(vrad_array);
      Py_XDECREF(sig_array);
      return NULL;
    }

    // Get array length
    int N = (int)PyArray_DIM(t_array, 0);

    /* Get pointers to the data as C-types. */
    double *t    = (double*)PyArray_DATA(t_array);
    double *vrad = (double*)PyArray_DATA(vrad_array);
    double *sig  = (double*)PyArray_DATA(sig_array);

    /**********************************************************************************/
    /* Start the complete clean algorithm, following the steps in Roberts et al. 1987 */
    /**********************************************************************************/
    double df, gain, SNlimit, span, ofac;
    int Nf, Nb, i, Nftot, Nwtot, iter;
    double *R_r, *R_i, *C_r, *C_i, *W_r, *W_i, *B_r, *S_r, *S_i;

    /* Define stopping criteria, frequencies, ... etc. */
    gain = 0.1; 	// value between 0.1 and 1
    SNlimit = 1.0; 	// procedure will continue until the last peak is less than SNlimit times the noise
    ofac = 4;     	// oversampling factor
    span = t[N-1] - t[0];	// Timespan
    df = 1.0 / (span*ofac);	// Frequency step
    Nf = floor((span*ofac)/plow);	// Amount of frequencies

    /* Define the total amount of frequencies - negative, positive + 0 */
    Nftot = 2*Nf + 1;
    Nwtot = 4*Nf + 1;
    /* Allocate memory for the spectra */
    R_r = (double *) malloc(Nftot * sizeof(double));
    R_i = (double *) malloc(Nftot * sizeof(double));
    C_r = (double *) malloc(Nftot * sizeof(double));
    C_i = (double *) malloc(Nftot * sizeof(double));
    W_r = (double *) malloc(Nwtot * sizeof(double));
    W_i = (double *) malloc(Nwtot * sizeof(double));
    S_r = (double *) malloc(Nftot * sizeof(double));
    S_i = (double *) malloc(Nftot * sizeof(double));

    /* Compute the dirty spectrum */
    FDFT(t, vrad, sig, N, df, Nf, R_r, R_i);

    /* Compute the window function */
    FDFT_Wfunc(t, N, df, 2*Nf, W_r, W_i);

    /* Do clean */
    doclean(R_r, R_i, W_r, W_i, Nf, C_r, C_i, gain, SNlimit, df, &iter);  // Steps 1,2,3 until we reach convergence
    GetCleanBeam(W_r, W_i, Nf, &B_r, &Nb); // Step 4
    Convolve_CleanBeam(C_r, C_i, R_r, R_i, S_r, S_i, Nf, Nb, B_r); // Step 5

    /* Compute the dirty spectrum again! */
    FDFT(t, vrad, sig, N, df, Nf, R_r, R_i);

    // Create python arrays

    npy_intp dim1[1], dim2[1];
    dim1[0] = 2*Nf+1;
    dim2[0] = 4*Nf+1;

    PyArrayObject *p_Rr;
    PyArrayObject *p_Ri;
    PyArrayObject *p_Wr;
    PyArrayObject *p_Wi;
    PyArrayObject *p_Sr;
    PyArrayObject *p_Si;

    p_Rr = (PyArrayObject *) PyArray_SimpleNew(1, dim1, PyArray_DOUBLE);
    p_Ri = (PyArrayObject *) PyArray_SimpleNew(1, dim1, PyArray_DOUBLE);
    p_Wr = (PyArrayObject *) PyArray_SimpleNew(1, dim2, PyArray_DOUBLE);
    p_Wi = (PyArrayObject *) PyArray_SimpleNew(1, dim2, PyArray_DOUBLE);
    p_Sr = (PyArrayObject *) PyArray_SimpleNew(1, dim1, PyArray_DOUBLE);
    p_Si = (PyArrayObject *) PyArray_SimpleNew(1, dim1, PyArray_DOUBLE);

    for (i = 0; i < 2*Nf+1; i++)
    {
        *(double *) (p_Rr->data + i*p_Rr->strides[0]) = R_r[i];
        *(double *) (p_Ri->data + i*p_Ri->strides[0]) = R_i[i];
        *(double *) (p_Sr->data + i*p_Sr->strides[0]) = S_r[i];
        *(double *) (p_Si->data + i*p_Si->strides[0]) = S_i[i];
    }
    
    for (i = 0; i < 4*Nf+1; i++)
    {
        *(double *) (p_Wr->data + i*p_Wr->strides[0]) = W_r[i];
        *(double *) (p_Wi->data + i*p_Wi->strides[0]) = W_i[i];
    }

    /* Clean up. */
    Py_DECREF(t_array);
    Py_DECREF(vrad_array);
    Py_DECREF(sig_array);
    free(R_r);
    free(R_i);
    free(C_r);
    free(C_i);
    free(W_r);
    free(W_i);
    free(S_r);
    free(S_i); 

    /* Build the output tuple */
    return Py_BuildValue("diiNNNNNN", df, Nf, iter, p_Rr, p_Ri, p_Wr, p_Wi, p_Sr, p_Si);

}

// docstring
static char clean_docs[] = 
"Run the CLEAN algorithm\n\n"
"Parameters\n"
"----------\n"
"time : array_like\n"
"       Array containing times of observations\n"
"rv : array_like\n"
"     Array containing the radial velocity values"
"err : array_like\n"
"      Array with associated uncertainties\n"
"plow : scalar\n"
"       Lowest period to sample.\n"
"\n"
"Returns\n"
"-------\n"
"df : scalar\n"
"     Frequency step\n"
"Nf : scalar\n"
"     Number of positive frequencies in the spectrum\n"
"iter : int\n"
"       Number of iterations performed \n"
"D_r, D_i : ndarrays\n"
"           Real and imaginary parts of the dirty spectrum\n"
"W_r, W_i : ndarrays\n"
"           Real and imaginary parts of the window function\n"
"C_r, C_i : ndarrays\n"
"           Real and imaginary parts of the cleaned spectrum\n"
"\n"
"Notes\n"
"-----\n"
"Algorithm from Roberts et al. 1987\n";


static PyMethodDef clean_methods[] = {
        {"clean", clean, METH_VARARGS, clean_docs},
        {NULL, NULL, 0, NULL}
};


void initperiodogram_CLEAN(void)
{    
        (void) Py_InitModule("periodogram_CLEAN", clean_methods);
        import_array();
            
}

#undef DFTTINY
#undef BEAMTOL

