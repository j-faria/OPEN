Tutorial
========

When you start OPEN you're presented with the shell prompt::


	$ python open.py
	Welcome to OPEN v0.0.3
	Created by João Faria | joao.faria@astro.up.pt

	Type "listcommands" for a list of the available commands
	For help on a specific command, type "command -h" or see http://j-faria.github.io/OPEN
	OPEN [1]: 

You can follow the instructions and type ``listcommands`` to see the available commands (this only lists the tested ones, everything else is still experimental). Whenever you type ``command -h`` you will see the usage pattern and options specific to that command. Each command's documentation is also available :doc:`here <commands>`

.. tip:: 
   Tab-completion will work on most OPEN commands, although it will automatically insert the leading percent sign (which is fine).

Let's start by reading radial velocity measurements from a sample file. 