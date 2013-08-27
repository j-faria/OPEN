# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
"""
This module defines the commands that are used as magics in OPEN. 
"""
from __future__ import absolute_import
from IPython.core.magic import (Magics, magics_class, line_magic, 
	                            needs_local_scope)
from IPython.core.magic_arguments import argument
from IPython.utils.io import ask_yes_no

from docopt import docopt
from classes import rvSeries
from logger import clogger, logging

################################################################################
################################################################################
## command usage patterns
read_usage = \
"""
Usage:
	read <file>...
	read <file>... -d
Options:
    -d 
"""

plot_usage = \
"""
Usage:
	plot [obs]
	plot [obs] -n SYSTEM
Options:
	-n SYSTEM 	  Specify name of system (else use default)
"""


# this is the most awesome function
def do_awesome(s):
    print 'I am a function and I do someting awesome!!!'
    print 'I also print this:', s

    
    
# These are additional magics that are exposed (only?) in embedded shells.
@magics_class
class EmbeddedMagics(Magics):

    @line_magic
    def kill_embedded(self, parameter_s=''):
        """%kill_embedded : deactivate for good the current embedded IPython.

        This function (after asking for confirmation) sets an internal flag so
        that an embedded IPython will never activate again.  This is useful to
        permanently disable a shell that is being called inside a loop: once
        you've figured out what you needed from it, you may then kill it and
        the program will then continue to run without the interactive shell
        interfering again.
        """

        kill = ask_yes_no("Are you sure you want to kill this embedded instance "
                         "(y/n)? [y/N] ",'n')
        if kill:
            self.shell.embedded_active = False
            print ("This embedded IPython will not reactivate anymore "
                   "once you exit.")

    @line_magic
    def awesome(self, parameter_s=''):
        do_awesome(parameter_s)

    @needs_local_scope
    @line_magic
    def read(self, parameter_s='', local_ns=None):
    	args = parse_arg_string('read', parameter_s)
    	filenames = args['<file>']

    	if local_ns.has_key('default') and not args['-d']:
    		return rvSeries(*filenames)
    	else:
    		local_ns['default'] = rvSeries(*filenames)

    @needs_local_scope
    @line_magic
    def plot(self, parameter_s='', local_ns=None):
    	args = parse_arg_string('plot', parameter_s)
    	print args

    	# use default system or user defined
    	try:
    		if local_ns.has_key('default') and not args['-n']:
    			system = local_ns['default']
    		else:
    			system_name = args['-n']
    			system = local_ns[system_name]
    	except KeyError:
    		from shell_colors import red
    		msg = red('ERROR: ') + 'Set a default system or provide a system '+\
    		                       'name with the -n option'
    		clogger.fatal(msg)
    		return

    	
    	if args['obs']:
    		system.do_plot_obs()


def parse_arg_string(command, arg_string):
	splitted = str(arg_string).split()

	if command is 'read':
		doc = read_usage
		args = docopt(doc, splitted)

	if command is 'plot':
		doc = plot_usage
		args = docopt(doc, splitted)

	return args
