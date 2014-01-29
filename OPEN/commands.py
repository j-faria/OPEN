# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
"""
This module defines the commands that are used as magics in OPEN. 
"""
# standard library imports
import glob
from itertools import chain

# IPython imports
from IPython.core.magic import (Magics, magics_class, line_magic, 
                                needs_local_scope)
from IPython.core.magic_arguments import argument
from IPython.utils.io import ask_yes_no

# intra-package imports
from .docopt import docopt, DocoptExit
from .classes import rvSeries
from .utils import stdout_write, selectable_plot, write_yorbit_macro
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
    read <file>... [-d] [--skip=<sn>] [--format=<form>] [-v]
    read -h | --help
Options:
    -d                  Set this as default system.
    -v --verbose        Verbose output about data just read.
    --skip=<sn>         How many header lines to skip [default: 0].
    --format=<form>     File format. One of 'drs35', 'drs34', 'coralie' or 'basic'.
    -h --help           Show this help message.
"""

plot_usage = \
"""
Usage:
    plot obs
    plot (fwhm | rhk | s | bis | contrast)
    plot -n SYSTEM
    plot -h | --help
Options:
    -n SYSTEM   Specify name of system (else use default)
    -h --help   Show this help message
"""

per_usage = \
"""
Usage:
    per obs
    per (bis | fwhm)
    per -n SYSTEM
    per (obs | bis | fwhm | resid) [--gls|--bayes|--fast] [-v] [--force] [--hifac=<hf>] [--fap]
    per -h | --help
Options:
    -n SYSTEM     Specify name of system (else use default)
    --gls         Calculate the Generalized Lomb-Scargle periodogram (default)
    --bayes       Calculate the Bayesian periodogram
    --fast        Calculate the Lomb-Scargle periodogram with fast algorithm
    --force       Force recalculation
    --hifac=<hf>  hifac * Nyquist is lowest frequency used [default: 40]
    --fap         Plot false alarm probabilities
    -v --verbose  Verbose statistical output 
    -h --help   Show this help message
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

dawfab_usage = \
"""
Usage:
    dawfab 
    dawfab -n SYSTEM
    dawfab -h | --help
Options:
    -n SYSTEM     Specify name of system (else use default)
    -h --help     Show this help message
"""

fit_usage = \
"""
Usage:
    fit [-v]
Options:
    -v --verbose  Verbose statistical output 
"""


correlate_usage = \
"""
Usage:
    correlate <var1> <var2> [-v] 
Options:
    -v --verbose  Verbose statistical output 
"""


restrict_usage = \
"""
Usage:
    restrict [(err <maxerr>)]
    restrict [(jd <minjd> <maxjd>)]
    restrict [(year <yr>)]
    restrict [(years <yr1> <yr2>)]
    restrict --gui
Options:
    --gui         Restrict data using a graphical interface (experimental)
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

        try:
            args = parse_arg_string('read', parameter_s)
        except DocoptExit:
            print read_usage.lstrip()
            return
        except SystemExit:
            return

        # take care of glob expansions
        globs = [glob.glob(f) for f in args['<file>']]
        filenames = list(chain.from_iterable(globs)) # some magic...

        # if 'default' system is already set, return the rvSeries class
        # this is useful when working with various systems simultaneously so 
        # that we can do, e.g., HDXXXX = %read file1 file2
        if local_ns.has_key('default') and not args['-d']:
            return rvSeries(*filenames, skip=args['--skip'])
        else:
            local_ns['default'] = rvSeries(*filenames, verbose=args['--verbose'],
                                                       skip=args['--skip'],
                                                       format=args['--format'])

    @needs_local_scope
    @line_magic
    def plot(self, parameter_s='', local_ns=None):
        """ Plot various quantities. 
        Type 'plot -h' for more help """

        try:
            args = parse_arg_string('plot', parameter_s)
        except DocoptExit:
            print plot_usage.lstrip()
            return
        except SystemExit:
            return
        # print args

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
        
        # plot the observed radial velocities
        if args['obs']:
            system.do_plot_obs()

        # plot other quantities
        extras_available = ['fwhm', 'contrast', 'bis_span', 'noise', 
                          's_mw', 'sig_s', 'rhk', 'sig_rhk', 
                          'sn_CaII', 'sn10', 'sn50', 'sn60']
        extras_mapping = ['fwhm', 'contrast', 'bis', 'noise',
                          's', 's', 'rhk', 'rhk',
                          'sn', 'sn', 'sn', 'sn']
        for i, e in enumerate(extras_available):
            try:
                if args[extras_mapping[i]]: 
                    system.do_plot_extras(e)
                    return
            except KeyError:
                pass

    @needs_local_scope
    @line_magic
    def per(self, parameter_s='', local_ns=None):
        """ Calculate periodograms of various quantities. 
        Type 'per -h' for more help. """

        try:
            args = parse_arg_string('per', parameter_s)
        except DocoptExit:
            print per_usage.lstrip()
            return
        except SystemExit:
            return
        #print args
        
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
        hf = float(args.pop('--hifac'))

        # which periodogram should be calculated?
        per_fcn = None
        if args['--bayes']: 
            per_fcn = periodograms.bls
            name = 'Bayesian Lomb-Scargle'
        if args['--fast']: 
            per_fcn = periodograms.LombScargle
            name = 'Lomb Scargle'
        if args['--gls']: 
            per_fcn = periodograms.gls
            name ='Generalized Lomb-Scargle'
        # this is the default if user did not specify arguments
        if per_fcn is None: 
            per_fcn = periodograms.gls
            name ='Generalized Lomb-Scargle'

        if args['obs']: # periodogram of the observed RVs
            try: 
                # are we forced to recalculate it?
                if args['--force']: raise AttributeError
                # it was calculated already?
                system.per 
                # the same periodogram?
                if system.per.name != name:
                    raise AttributeError
                # system.per._output(verbose=verb)  # not ready
                system.per._plot(doFAP=args['--fap'])
            except AttributeError:
                system.per = per_fcn(system, hifac=hf)
                # system.per._output(verbose=verb)  # not ready
                system.per._plot(doFAP=args['--fap'])
                print system.per.name

        if args['bis']: # periodogram of the CCF's Bisector Inverse Slope
            system.bis_per = per_fcn(system, quantity='bis')
            system.bis_per._plot()

        if args['fwhm']: # periodogram of the CCF's fwhm
            system.fwhm_per = per_fcn(system, quantity='fwhm')
            system.fwhm_per._plot()

        if args['resid']: # periodogram of the residuals of the current fit
            system.resid_per = per_fcn(system, quantity='resid')
            system.resid_per._plot()

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


    @needs_local_scope
    @line_magic
    def dawfab(self, parameter_s='', local_ns=None):
        """ Run the Dawson Fabrycky algorithm to search for aliases.
        Type 'dawfab -h' for more help. """

        args = parse_arg_string('dawfab', parameter_s)
        if args == 1: return
        #print args
        
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
        
        core.do_Dawson_Fabrycky(system)
        return


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
        #print args

        verb = True if args['--verbose'] else False

        if local_ns.has_key('default'):
            system = local_ns['default']
            result = core.do_fit(system, verb)
            if result is not None:
                system.model['drift'] = result
                if verb: print 'Coeff:', result
                system.do_plot_drift()

        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        system.do_plot_fit()

    @needs_local_scope
    @line_magic
    def correlate(self, parameter_s='', local_ns=None):
        from shell_colors import red
        args = parse_arg_string('correlate', parameter_s)
        if args == 1: return
        #print args

        verb = True if args['--verbose'] else False

        if local_ns.has_key('default'):
            system = local_ns['default']
            var1 = args['<var1>']
            var2 = args['<var2>']
            result = core.do_correlate(system, vars=(var1, var2), verbose=verb)
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

    @needs_local_scope
    @line_magic
    def gen(self, parameter_s='', local_ns=None):
        """ Run the genetic algorithm minimization - stub """
        from shell_colors import red
        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        core.do_genetic(system)
        system.do_plot_fit()

    @needs_local_scope
    @line_magic
    def genyorbit(self, parameter_s='', local_ns=None):
        """ Run the genetic algorithm minimization - stub """
        from shell_colors import red
        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return
            
        write_yorbit_macro(system)
        # core.do_genetic(system)
        # system.do_plot_fit()

    @needs_local_scope
    @line_magic
    def nest(self, parameter_s='', local_ns=None):
        """ Start the MultiNest analysis and handle data interaction and IO """
        from shell_colors import red

        if local_ns.has_key('default'):
            system = local_ns['default']
            core.do_multinest(system)
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

        if args == 1: return
        print args

        # use default system or user defined
        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return


        if args['err']: 
            try:
                maxerr = float(args['<maxerr>'])
            except ValueError:
                msg = red('ERROR: ') + 'maxerr shoud be a number!'
                clogger.fatal(msg)
                return
            core.do_restrict(system, 'error', maxerr)

        if args['jd']:
            try:
                maxjd = int(args['<maxjd>'])
                minjd = int(args['<minjd>'])
            except ValueError:
                msg = red('ERROR: ') + 'minjd and maxjd shoud be numbers!'
                clogger.fatal(msg)
                return
            core.do_restrict(system, 'date', minjd, maxjd)

        if args['year']:
            try:
                yr = int(args['<yr>'])
            except ValueError:
                msg = red('ERROR: ') + 'yr shoud be a number!'
                clogger.fatal(msg)
                return
            core.do_restrict(system, 'year', yr)

        if args['years']:
            try:
                yr1 = int(args['<yr1>'])
                yr2 = int(args['<yr2>'])
            except ValueError:
                msg = red('ERROR: ') + 'yr1 and yr2 shoud be numbers!'
                clogger.fatal(msg)
                return
            core.do_restrict(system, 'years', yr1, yr2)

        if args['--gui']:
            selectable_plot([1,2,3], [4, 16, 32], 'ro')





def parse_arg_string(command, arg_string):
    """ Parse arguments for each of the commands. """
    # docopt does the heavy-lifting parsing, we just split the argument string
    # and catch the exceptions raised by -h or --help

    splitted = str(arg_string).split()

    if command is 'read':
        args = docopt(read_usage, splitted)

    if command is 'plot':
        args = docopt(plot_usage, splitted)

    if command is 'per':
        args = docopt(per_usage, splitted)

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

    if command is 'correlate':
        try:
            args = docopt(correlate_usage, splitted)
        except SystemExit:
            return 1

    if command is 'restrict':
        if arg_string == '': 
            return DocoptExit # restrict needs arguments
        try:
            args = docopt(restrict_usage, splitted)
        except (SystemExit, DocoptExit) as e:
            return 1

    if command is 'wf':
        try:
            args = docopt(wf_usage, splitted)
        except SystemExit:
            return 1

    if command is 'dawfab':
        try:
            args = docopt(dawfab_usage, splitted)
        except SystemExit:
            return 1

    return args
