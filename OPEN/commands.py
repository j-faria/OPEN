# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
"""
This module defines the commands that are used as magics in OPEN. 
"""
# IPython imports
from IPython.core.magic import (Magics, magics_class, line_magic, 
                                needs_local_scope)
from IPython.core.magic_arguments import argument
from IPython.utils.io import ask_yes_no

# intra-package imports
from .docopt import docopt, DocoptExit
from .classes import rvSeries
from .utils import stdout_write
from .logger import clogger, logging
import core
import periodograms

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
    -n SYSTEM       Specify name of system (else use default)
"""

per_usage = \
"""
Usage:
    per [obs]
    per [obs] -n SYSTEM
    per [obs] -v
    per -h | --help
Options:
    -n SYSTEM     Specify name of system (else use default)
    -v --verbose  Verbose statistical output 
    -h --help     Show this help message
"""


wf_usage = \
"""
Usage:
    wf 
    wf -n SYSTEM
    wf --dials
    wf -h | --help
Options:
    -n SYSTEM     Specify name of system (else use default)
    --dials       Plot phase "dials" in largest (3) peaks
    -h --help     Show this help message
"""


fit_usage = \
"""
Usage:
    fit [-v]
Options:
    -v --verbose  Verbose statistical output 
"""


restrict_usage = \
"""
Usage:
    restrict [(err <maxerr>)]
Options:
    -h --help     Show this help message
"""

    
command_list = \
"""
 read       Read RV files. 
 plot       Plot various quantities.
 per        Calculate periodograms.
 mod        Define the model that will be adjusted to the data.
 fit        ...
 restrict   Select data based on date, SNR or RV accuracy.
"""


    
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


    @needs_local_scope
    @line_magic
    def read(self, parameter_s='', local_ns=None):
        """ Read files with RV measurements.
        Type 'read -h' for more help """

        args = parse_arg_string('read', parameter_s)
        filenames = args['<file>']

        # if 'default' system is already set, return the rvSeries class
        # this is useful when working with various systems simultaneously so 
        # that we can do, e.g., HDXXXX = %read file1 file2
        if local_ns.has_key('default') and not args['-d']:
            return rvSeries(*filenames)
        else:
            local_ns['default'] = rvSeries(*filenames)

    @needs_local_scope
    @line_magic
    def plot(self, parameter_s='', local_ns=None):
        """ Plot various quantities. 
        Type 'plot -h' for more help """
        args = parse_arg_string('plot', parameter_s)
        if args == 1: return
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

    @needs_local_scope
    @line_magic
    def per(self, parameter_s='', local_ns=None):
        """ Calculate periodograms of various quantities. 
        Type 'per -h' for more help. """

        args = parse_arg_string('per', parameter_s)
        if args == 1: return
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
        
        verb = True if args['--verbose'] else False
        if args['obs']:
            try: 
                system.per
                system.per._output(verbose=verb)
                system.per._plot()
            except AttributeError:
                system.per = periodograms.gls(system, hifac=5)
                system.per._output(verbose=verb)
                system.per._plot()

    @needs_local_scope
    @line_magic
    def wf(self, parameter_s='', local_ns=None):
        """ Calculate the spectral window function of the observations. 
        Type 'wf -h' for more help. """

        args = parse_arg_string('wf', parameter_s)
        if args == 1: return
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
        
        try: 
            system.per
        except AttributeError:
            from shell_colors import green
            clogger.debug('Calculating periodogram to get frequencies')
            stdout_write('Calculating periodogram to get frequencies...')
            system.per = periodograms.gls(system, hifac=5)
            print green(' done')
        
        try: 
            system.wf._plot()
        except AttributeError:
            system.wf = periodograms.SpectralWindow(system.per.freq, system.time)

    @line_magic
    def listcommands(self, parameter_s=''):
        """ List available commands """
        print command_list

    @needs_local_scope
    @line_magic
    def mod(self, parameter_s='', local_ns=None):
        """ Define the type of model that will be adjusted to the data.
        Type 'mod -h' for more help
        """
        from shell_colors import yellow, blue, red
        args = parse_arg_string('mod', parameter_s)
        if args == 1:  # called without arguments, show how it's done
            msg = yellow('Usage: ') + 'mod [k<n>] [d<n>]\n' + \
                  'Options: k<n>    Number of keplerian signals\n' + \
                  '         d<n>    Degree of polynomial drift'
            clogger.fatal(msg)
            return

        if local_ns.has_key('default'):
            system = local_ns['default']
            system.model = {}
            system.model['k'] = k = int(args[0][1])
            system.model['d'] = d = int(args[1][1])
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        # this should be logged?
        print blue('Current model:'), k, 'kep,', d, 'drifts'

        # ... do someting with this ...

    @needs_local_scope
    @line_magic
    def fit(self, parameter_s='', local_ns=None):
        from shell_colors import red
        args = parse_arg_string('fit', parameter_s)
        if args == 1: return
        print args

        verb = True if args['--verbose'] else False

        if local_ns.has_key('default'):
            system = local_ns['default']
            result = core.do_fit(system, verb)
            if result is not None:
                system.model['drift'] = result
                if verb: print 'Coeff:', result
                system.do_plot_drift()
            else:
                return
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

    @needs_local_scope
    @line_magic
    def gen(self, parameter_s='', local_ns=None):
        from shell_colors import red

        if local_ns.has_key('default'):
            system = local_ns['default']
            core.do_genetic(system)
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

    @needs_local_scope
    @line_magic
    def restrict(self, parameter_s='', local_ns=None):
        """ Select data based on date, SNR or radial velocity accuracy.
        Type 'restrict -h' for more help
        """
        from shell_colors import yellow, blue, red
        args = parse_arg_string('restrict', parameter_s)
        if args == DocoptExit:
            msg = yellow('Warning: ') + "I'm not doing anything. See restrict -h"
            clogger.fatal(msg)
            return

        print args

        # use default system or user defined
        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            print 'No default!'  # handle this!
            return


        if args['err']: 
            try:
                maxerr = int(args['<maxerr>'])
            except ValueError:
                msg = red('ERROR: ') + 'maxerr shoud be a number!'
                clogger.fatal(msg)
                return
            core.do_restrict(system, 'error', maxerr)





def parse_arg_string(command, arg_string):
    """ Parse arguments for each of the commands. """
    # docopt does the heavy-lifting parsing, we just split the argument string
    # and catch the exceptions raised by -h or --help

    splitted = str(arg_string).split()

    if command is 'read':
        doc = read_usage
        args = docopt(doc, splitted)

    if command is 'plot':
        try:
            args = docopt(plot_usage, splitted)
        except SystemExit:
            return 1

    if command is 'per':
        try:
            args = docopt(per_usage, splitted)
        except SystemExit:
            return 1

    # this is a little different
    if command is 'mod':
        import re
        if arg_string == '': 
            return 1 # mod needs arguments
        if arg_string in ('-h', '--help'):
            return 1 # explain what to do

        k = re.compile("k[0-9]").findall(arg_string)
        if k == []: # if just drifts
            k = ['k0']
        d = re.compile("d[0-9]").findall(arg_string)
        if d == []: # if just keplerians
            d = ['d0']
        args = k+d

    if command is 'fit':
        try:
            args = docopt(fit_usage, splitted)
        except SystemExit:
            return 1

    if command is 'restrict':
        if arg_string == '': 
            return DocoptExit # restrict needs arguments
        try:
            args = docopt(restrict_usage, splitted)
        except (SystemExit, DocoptExit) as e:
            return e

    if command is 'wf':
        try:
            args = docopt(wf_usage, splitted)
        except SystemExit:
            return 1

    return args
