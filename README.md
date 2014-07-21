OPEN
====

Tools for exoplanet detection.

- A description of your project
- Links to the project's ReadTheDocs page
- A TravisCI button showing the state of the build
- "Quickstart" documentation (how to quickly install and use your project)
- A list of non-Python dependencies (if any) and how to install them
                    

Quickstart
----------

OPEN depends on the following software:

   * [Numpy](http://www.numpy.org/)
   * [matplotlib](http://matplotlib.org/)
   * [IPython](http://ipython.org/)
 
   * Other non-Python packages

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
