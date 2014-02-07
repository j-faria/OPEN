OPEN
====

Tools for exoplanet detection.

- A description of your project
- Links to the project's ReadTheDocs page
- A TravisCI button showing the state of the build
- "Quickstart" documentation (how to quickly install and use your project)
- A list of non-Python dependencies (if any) and how to install them
                    

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
