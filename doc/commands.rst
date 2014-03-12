.. _commands-page:

Commands
========

At its core, OPEN is just an IPython shell enhanced with a few custom magic functions. You might be familiar with magic functions, things like ``%cd``, ``%run`` or ``%timeit``. These can often be called without the leading '%' from within IPython and offer a convenient way to execute common functions which have always the same type of parameters (or no parameters at all). Thus, by creating new magic functions we are basically adding custom commands to the shell. This is in addition to everything that is built-in to both Python and IPython so that all batteries are included; all the power of Python is there.


Currently, these are the commands available in `OPEN` 

read
~~~~
Read files with RV measurements.

.. program-output:: /home/joao/anaconda/bin/python -c "import os; os.chdir('..'); from OPEN.commands import read_usage; print read_usage"


plot
~~~~
Plot various quantities.

.. program-output:: /home/joao/anaconda/bin/python -c "import os; os.chdir('..'); from OPEN.commands import plot_usage; print plot_usage"


per
~~~
Calculate periodograms.

.. program-output:: /home/joao/anaconda/bin/python -c "import os; os.chdir('..'); from OPEN.commands import per_usage; print per_usage"


mod
~~~
Define the model that will be adjusted to the data. ::

	Usage:
	    mod [k<n>] [d<n>]
	Options: 
	    k<n>    Number of keplerian signals
	    d<n>    Degree of polynomial drift


restrict
~~~~~~~~
Select data based on date, SNR or RV accuracy.

.. program-output:: /home/joao/anaconda/bin/python -c "import os; os.chdir('..'); from OPEN.commands import restrict_usage; print restrict_usage"


killall
~~~~~~~
Close all plot windows.

