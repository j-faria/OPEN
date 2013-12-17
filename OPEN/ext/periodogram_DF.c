#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#define TWOPID 6.2831853071795865 // 2 * pi

// Wordt aangeroepen door het python programma
static PyObject *periodogram_DF(PyObject *self, PyObject *args)
{

    /* constants */
    PyArrayObject *x;           // time
    PyArrayObject *y;           // RV
    PyArrayObject *sig = NULL;  // errors

    double ofac = 1;            // Oversamplingfactor
    double plow = 0.5;          // Lowest period to sample
    double median_step;         // Median step in time
    double xmax;                // Maximum of x
    double xmin;                // Minimum of x
    double xdif;                // Range in x
    double fNy;                 // Nyquist frequency

    double temp;                // Temporary value

    int n;                      // Size of array
    int i;
    int j;

    int nout;                   // Amount of frequencies


    /* Everything after | is optional */
    if ( !PyArg_ParseTuple(args, "O!O!|O!dd",
            &PyArray_Type, &x,
            &PyArray_Type, &y,
            &PyArray_Type, &sig,
            &ofac,&plow)
       )
    {
        printf("potjepotje");
        return NULL;
    }

    // Set array length
    n = x->dimensions[0];

    // Create array of weights
    double weight[n];

    // Sig is optional
    if (sig == NULL)
    {
        for (i = 0; i < n; i++)
        {
            weight[i] = 1 / n;
        }
    } else {
        double W = 0.;
        double k;      // double from sig array 
        
        for (i = 0; i < n; i++)
        {
            k = *(double *) (sig->data + i*sig->strides[0]);
            W += 1 / (k * k);
        }

        for (i = 0; i < n; i++)
        {
            k = *(double *) (sig->data + i*sig->strides[0]);
            weight[i] = 1 / W * 1 / (k * k);
        }
    }
    

    // Calculate time differences
    double xsteps[n-1];

    for (i = 0; i < n-1; i++)
    {
        xsteps[i] = *(double *) (x->data + (i+1)*x->strides[0]) - *(double *) (x->data + i*x->strides[0]);
    }

    // Calculate median
    for (i = 0; i < n-2; i++)
    {
        for (j = i+1; j < n-1; j++)
        {
            if (xsteps[j]<xsteps[i])
            {
                temp = xsteps[i];
                xsteps[i]=xsteps[j];
                xsteps[j]=temp;
            }
        }
    }

    if ( (n-1)%2 == 0 )
    {
        median_step = ((double) ((xsteps[(n-1)/2]+xsteps[(n-1)/2-1])/2.0));
    }
    else
    {
        median_step = ((double) (xsteps[(n-1)/2]));
    }

//    if (median_step < 0.4)
//    {
//        median_step=2.0;
//    }

    // Nyquist frequency

    fNy = 1. / (2. * median_step);

    // Calculate the minimum and maximum of x
    xmin = *(double *) (x->data);
    xmax = *(double *) (x->data);
    
    for (i = 0; i < n; i++)
    {
        temp = *(double *) (x->data + i*x->strides[0]);
        if (temp > xmax)
        {
            xmax = temp;
        }

        if (temp < xmin)
        {
            xmin = temp;
        }
    }

    xdif = xmax - xmin;


    // Calculate the amount of frequencies
//    nout = (int) (xdif * ofac / (2. * median_step) - 1);     // ending with Nyquist
    nout = (int) (xdif * ofac / plow);      // ending with the lowest period requested

    double freqnow;
    double arg;

    double wpr[n];
    double wpi[n];
    double wr[n];
    double wi[n];

    double px[nout];
    double py[nout];

    double a_cos[nout];
    double b_sin[nout];
    double c_cte[nout];
    double phi[nout];

    double s;
    double c;
    double sums;
    double sumc;
    double sumss;
    double sumcc;
    double sumcs;
    double sumy;
    double yy;
    double ydif;
    double ys;
    double yc;
    double cc;
    double ss;
    double cs;
    double d;
    double wtemp;
    



    /* Get the powers for each frequency and the maximum */

        freqnow=1.0/(xdif*ofac); /* start frequency = step in frequencies */


        for (j=0;j<n;j++) {        /* Step through time */
                /* initialise the values for the trigonometric recurrences */
                temp = *(double *) (x->data + j*x->strides[0]);
                arg=TWOPID*((temp)*freqnow); /* arg = 2*pi*t*f = omega*t */
                wpr[j] = -2.0*sin(0.5*arg)*sin(0.5*arg);
                wpi[j]=sin(arg); 
                wr[j]=cos(arg); 
                wi[j]=wpi[j];   
        }

        for (i = 0; i < nout; i++) {        /* Step through frequency */
                /* Store frequency */
                px[i]=freqnow;
                /* Follow Zechmeister 2009 */
                sums=sumc=sumss=sumcc=sumcs=sumy=ydif=yy=ys=yc=0.0;
                for (j=0;j<n;j++) {
                        temp = *(double *) (y->data + j*y->strides[0]);
                        c=wr[j]; /* c=cos(omega*t) */
                        s=wi[j]; /* s=sin(omega*t) */
                        sums += weight[j]*s;    /* S */
                        sumc += weight[j]*c;    /* C */
                        sumss += weight[j]*s*s; /* hat SS */
                        sumcc += weight[j]*c*c; /* hat CC */
                        sumcs += weight[j]*s*c; /* hat CS */
                        sumy += weight[j]*temp;  /* Y = y bar */
                        
                }    
                for (j=0;j<n;j++) {
                        temp = *(double *) (y->data + j*y->strides[0]);
                        c=wr[j]; /* c=cos(omega*t) */
                        s=wi[j]; /* s=sin(omega*t) */
                        ydif=temp-sumy;           
                        yy += weight[j]*ydif*ydif; /* YY */
                        ys += weight[j]*ydif*s; /* YS */
                        yc += weight[j]*ydif*c; /* YC */
                        /* Define the next cosines and sines with recurrence formula */
                        wtemp = wr[j];
                        wr[j]=(wtemp*wpr[j]-wi[j]*wpi[j])+wr[j];
                        wi[j]=(wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
                }

                cc = sumcc - sumc*sumc; /* CC */
                ss = sumss - sums*sums; /* SS */
                cs = sumcs - sumc*sums; /* CS */
                d = cc*ss - cs*cs; /* D */

                /* The power */      
                py[i]=1.0/(yy*d) * (ss*yc*yc + cc*ys*ys - 2.0*cs*ys*yc);

                /* Calculate the results */

                /* The result of the fit: a cos(omega t) + b sin(omega t) + c */
                a_cos[i] = (yc*ss-ys*cs)/d;
                b_sin[i] = (ys*cc-yc*cs)/d;
                c_cte[i] = sumy - a_cos[i]*sumc - b_sin[i]*sums;

                /* The phase */
                phi[i] = atan2(a_cos[i],b_sin[i]);

                /* Define the next frequency */
                freqnow += 1.0/(ofac*xdif);  
        }



    // Create python arrays
    int dim[1];
    dim[0] = nout;

    PyArrayObject *p_px;
    PyArrayObject *p_py;
    PyArrayObject *p_a_cos;
    PyArrayObject *p_b_sin;
    PyArrayObject *p_c_cte;
    PyArrayObject *p_phi;

    p_px = (PyArrayObject *) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
    p_py = (PyArrayObject *) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
    p_a_cos = (PyArrayObject *) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
    p_b_sin = (PyArrayObject *) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
    p_c_cte = (PyArrayObject *) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
    p_phi = (PyArrayObject *) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);

    for (i = 0; i < nout; i++)
    {
        *(double *) (p_px->data + i*p_px->strides[0]) = px[i];
        *(double *) (p_py->data + i*p_py->strides[0]) = py[i];
        *(double *) (p_a_cos->data + i*p_a_cos->strides[0]) = a_cos[i];
        *(double *) (p_b_sin->data + i*p_b_sin->strides[0]) = b_sin[i];
        *(double *) (p_c_cte->data + i*p_c_cte->strides[0]) = c_cte[i];
        *(double *) (p_phi->data + i*p_phi->strides[0]) = phi[i];

    }
    


    return Py_BuildValue("NNNNNNdd", p_px, p_py, p_a_cos, p_b_sin, p_c_cte, p_phi, fNy, xdif);

}

static PyMethodDef periodogram_DF_methods[] = {

        {"periodogram_DF", periodogram_DF, METH_VARARGS,
         "Compute the Lomb-Scargle periodogram + phases."},
         
        {NULL, NULL, 0, NULL}

};


void initperiodogram_DF(void)
{    
        (void) Py_InitModule("periodogram_DF", periodogram_DF_methods);
            
        import_array();
            
}



#undef TWOPID
