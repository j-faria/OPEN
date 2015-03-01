OPEN
====

Tools for exoplanet detection.


Dependencies
------------

OPEN depends on the following Python packages:

   * [Numpy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), [IPython](http://ipython.org/), [Scipy](http://scipy.org/)

Different versions of these packages should not be causing problems. Nevertheless, OPEN is being developed with Numpy 1.9.1, Matplotlib 1.4.2, IPython 2.2.0 and Scipy 0.14.0 and so these versions are the ones officially supported. Please [file an issue](https://github.com/j-faria/OPEN/issues) if you notice any problems with older versions.

OPEN also depends on a fairly recent version of the [GNU Compiler Collection](https://gcc.gnu.org/) (at least v4.7) to compile Fortran routines and Python extensions. An easy way to obtain all the compiler/library dependencies [Note: only for 64-bit systems] is to use the [MESA SDK](http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk) and specify the path to `gfortran` during installation (see below).


Quickstart
----------

The easiest way to download OPEN is by using `git`. Open a terminal and run

    git clone --recursive git@github.com:j-faria/OPEN.git

which will create a directory `OPEN` containing the code (and the MultiNest submodule). Currently there is no other way to obtain the MultiNest submodule. This is a [known issue](https://github.com/j-faria/OPEN/issues/6).

To install OPEN, go into the directory you just created and type

    $ [cd OPEN]
    $ python scons/scons.py

or, if you are using `gfortran` from the MESA SDK,

    $ python scons/scons.py --gfortran=/path/to/mesasdk/bin/gfortran

If everything goes smoothly you can now run OPEN:

    $ python open.py

Troubleshooting
---------------

If the build succeds without errors but when you run `python open.py` you get the error
  
    ImportError: /home/joao/OPEN/OPEN/ext/gp.so: undefined symbol: dpotrf_

you need to install LAPACK as a shared library (for example: `apt-get install liblapack-dev`) and re-build OPEN with the command

    $ python scons/scons.py --gp-with-lapack

