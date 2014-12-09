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




Implementation details
----------------------

At its core, OPEN is just an IPython shell enhanced with a few custom magic functions. 
You might be familiar with magic functions, things like ```%cd```, ```%run``` or ```%timeit```. These can
often be called without the leading '%' from within IPython and offer a convenient way 
to execute common functions which have always the same type of parameters (or no parameters
at all). 
Then, by creating new magic functions we are basically adding custom commands to the 
shell. This is **in addition** to everything that is built-in to both Python and IPython
so that all batteries are included; all the power of Python is there.



The basics
----------

When you start OPEN you're presented with the shell prompt

```
$ python open.py 
Welcome to OPEN v0.0.2
Created by João Faria | joao.faria@astro.up.pt

Type "listcommands" for a list of the available commands
For help on a specific command, type "command-h" or see http://j-faria.github.io/OPEN
OPEN [1]: 
```

You can follow the instructions and type ```listcommands``` to see the available commands
(this only lists the tested ones, everything else is still experimental).
Whenever you type ```command -h``` you will see the usage pattern and options specific to
that command. 

> > **TIP:** Tab-completion will work on most OPEN commands, although it will automatically
insert the leading percent sign (which is fine)

Let's start by reading radial velocity measurements from a sample file. 
TODO...
