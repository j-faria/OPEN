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
from os.path import expanduser
from itertools import chain

# IPython imports
from IPython.core.magic import (Magics, magics_class, line_magic, 
                                needs_local_scope)
from IPython.core.magic_arguments import argument

# other imports
from numpy import sqrt, mean, min, delete, take

# intra-package imports
from .docopt import docopt, DocoptExit
from .classes import rvSeries
from .utils import stdout_write, ask_yes_no, write_yorbit_macro
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
    read <file>... [-d] [--skip=<sn>] [-v] [--quiet] [--nomps]
    read -h | --help
Options:
    -d                  Set this as default system.
    -v --verbose        Verbose output about data just read.
    --quiet             Do not print any output.
    --skip=<sn>         How many header lines to skip [default: 0].
    --nomps             Do not convert data to m/s
    -h --help           Show this help message.
"""

saverdb_usage = \
"""
Usage:
    saverdb <file>
    saverdb -n SYSTEM
Options:
    -n SYSTEM   Specify name of system (else use default)
"""

plot_usage = \
"""
Usage:
    plot (obs | fwhm | rhk | s | bis | contrast | resid) [--save=filename]
    plot -n SYSTEM
    plot -h | --help
Options:
    -n SYSTEM        Specify name of system (else use default)
    --save=filename  Save figure as filename
    -h --help        Show this help message
"""

per_usage = \
"""
Usage:
    per 
    per [-n SYSTEM] (obs|bis|fwhm|rhk|contrast|resid) [-g|-m|-b|-l|-z] [-v] [-f] [--hifac=<hf>] [--ofac=<of>] [--fap] [--bfap] [--save=filename] [--noplot]
    per -h | --help
Options:
    -n SYSTEM        Specify name of system (else use default)
    -g --gls         Calculate the Generalized Lomb-Scargle periodogram (default)
    -m --bgls        Calculate the Bayesian Generalized Lomb-Scargle periodogram
    -b --bayes       Calculate the Bayesian LS periodogram
    -l --ls          Calculate the Lomb-Scargle periodogram with fast algorithm
    -z --hoef        Calculate the Hoeffding-test "periodogram" with Zucker's algorithm
    -f --force       Force recalculation
    --hifac=<hf>     hifac * Nyquist is lowest frequency used [default: 40]
    --ofac=<of>      Oversampling factor [default: 6]
    --fap            Plot false alarm probabilities
    --bfap           Plot false alarm probabilities calculated using bootstrap
    --save=filename  Save figure as filename
    --noplot         Don't plot the periodogram (just creates system.per* instance)
    -v --verbose     Verbose statistical output 
    -h --help        Show this help message
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
    correlate <var1> <var2> [-v] [-r]
Options:
    -v --verbose  Verbose statistical output 
    -r --remove   Remove linear dependence from RV
"""


de_usage = \
"""
Usage:
    de [--npop=<pop>] [--ngen=<gen>]
    de -h | --help
Options:
    --npop=<pop>  Number of individuals in population [default: 100]
    --ngen=<gen>  Number of generations to evolve [default: 250]
    -h --help     Show this help message
"""


demc_usage = \
"""
Usage:
    demc [<zipfile>]
    demc -n SYSTEM
Options:
    -n SYSTEM   Specify name of system (else use default)
"""


gen_usage = \
"""
Usage:
    gen [--npop=<pop>] [--ngen=<gen>]
    gen -h | --help
Options:
    --npop=<pop>  Number of individuals in population [default: 100]
    --ngen=<gen>  Number of generations to evolve [default: 250]
    -h --help     Show this help message
"""


rrot_usage = \
"""
Usage:
    remove_rotation [fwhm]
    remove_rotation [-n SYSTEM] [fwhm | rhk] [--prot=<p>] [--nrem=<nr>]
    remove_rotation -h | --help
Options:
    -n SYSTEM   Specify name of system (else use default)
    --prot=<p>
    --nrem=<nr>  Number of harmonics to remove (including Prot) [default: 1]
    -h --help    Show this help message
"""

addnoise_usage = \
"""
Usage:
    add_noise <number>
    add_noise -n SYSTEM <number>
    add_noise -h | --help
Options
    <number>    Add <number> m/s to the RV uncertainties.
    -n SYSTEM   Specify name of system (else use default)
    -h --help   Show this help message
"""

nest_usage = \
"""
Usage:
    nest 
    nest [options]
Options:
    -u            User sets the namelist file
    -r            Resume from a previous MultiNest run
    -v            Be verbose on output and plots
    --gp          Perform model selection within Gaussian Process framework
    --jitter      Include a jitter parameter (incompatible with --gp)
    --train=None  Train the GP on quantity before using it in the RVs
    --skip-mcmc   Skip the training MCMC: the user sets the appropriate namelist options
    --lin=None    Include linear dependence on quantity in the model
    --ncpu=<cpu>  Number of threads to use; by default use all available
    --noplot      Do not produce result plots
    --saveplot    Save all plots (does nothing if --noplot is given)
    --feed        Force feedback on the progress of the evidence calculation
    --MAPfeed     Force feedback on the current MAP parameters
    --maxp=mp     Maximum number of planets to include in automatic run [default: 3]
    --restart     Fully restart a previous automatic model selection run
    --nml=None    Specify the `full` path to the namelist file
    --startp=None Comma-separated list of planets to start from the beginning (overide -r)
"""

restrict_usage = \
"""
Usage:
    restrict [(err <maxerr>)]
    restrict [(jd <minjd> <maxjd>)]
    restrict [(year <yr>)]
    restrict [(years <yr1> <yr2>)]
    restrict --gui
    restrict --index=None [--noask]
Options:
    --gui         Restrict data using a graphical interface (experimental)
    --index=None  Remove specific data points, providing their indices [default:None]
    --noask       Do not confirm if removing observations
"""


rotation_usage = \
"""
Usage:
    rotation
    rotation -n SYSTEM
    rotation -h | --help
Options
    -n SYSTEM   Specify name of system (else use default)
    -h --help   Show this help message
"""

tp_mps_usage = \
"""
Usage:
    to_mps
    to_mps -n SYSTEM
    to_mps -h | --help
Options
    -n SYSTEM   Specify name of system (else use default)
    -h --help   Show this help message
"""

create_usage = \
"""
Usage:
    create
    create np(%d) [p(%f)] [e(%f)] [k(%f)] [N(%d)] [out(%s)] [sample(%s)]
    create --gui
Options:
    np(%d) [p(%f)] [e(%f)] [k(%f)] [N(%d)] [out(%s)] [sample(%s)]   Batch processing
    --gui    Create data using a graphical interface (experimental)
"""


    
command_list = \
"""
 read       Read RV files. 
 plot       Plot various quantities.
 per        Calculate periodograms.
 mod        Define the model that will be adjusted to the data.
 de         Fit the model using a Differential Evolution algorithm (somewhat experimental...)
 restrict   Select data based on date, SNR or RV accuracy.
 rotation   Calculate rotation period from activity-rotation relation.
 create     Generate synthetic RV data.
 killall    Close all plot windows
"""
    
# These are additional magics that are exposed (only?) in embedded shells.
@magics_class
class EmbeddedMagics(Magics):

    @line_magic
    def develop(self, parameter_s=''):
        # reload(classes)
        import reimport, os
        mod = reimport.modified()
        reimport.reimport(*mod)
        print 'Done re-importing'
        # reload(periodograms)
        # reload(core)
        os.system('python scons/scons.py --gfortran=/home/joao/Software/mesasdk/bin/gfortran')


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

        # take care of glob (and tilde) expansions
        files = args['<file>']

        # hack for metal-poor files
        if len(files) == 1 and files[0].startswith('HD'):
            files = ['/home/joao/phd/data/'+files[0]+'_harps_mean_corr.rdb']
        ##
        globs = [glob.glob(expanduser(f)) for f in files]
        filenames = list(chain.from_iterable(globs)) # some magic...

        # if 'default' system is already set, return the rvSeries class
        # this is useful when working with various systems simultaneously so 
        # that we can do, e.g., HDXXXX = %read file1 file2
        if local_ns.has_key('default') and not args['-d']:
            try:
                return rvSeries(*filenames, skip=args['--skip'], verbose=not args['--quiet'])
            except AttributeError:
                pass
        else:
            try:
                local_ns['default'] = rvSeries(*filenames, skip=args['--skip'], verbose=not args['--quiet'])
            except IOError:
                return

        default = local_ns['default']
        
        if args['--verbose'] and not args['--quiet']:
            default.stats()

        if (min(default.error) < 0.01 and not args['--nomps']):
            from shell_colors import blue
            mean_vrad = mean(default.vrad)

            if not args['--quiet']:
                # msg = blue('INFO: ') + 'Converting to m/s and subtracting mean value of %f' % mean_vrad
                msg = blue('INFO: ') + 'Converting to m/s'
                clogger.info(msg)

            default.vrad = (default.vrad - mean_vrad)*1e3 + mean_vrad
            default.error *= 1e3
            default.vrad_full = (default.vrad_full - mean(default.vrad_full))*1e3 + mean(default.vrad_full)
            default.error_full *= 1e3
            default.units = 'm/s'

    @needs_local_scope
    @line_magic
    def saverdb(self, parameter_s='', local_ns=None):
        """ Save current system's RV in a file """
        try:
            args = parse_arg_string('saverdb', parameter_s)
        except DocoptExit:
            print saverdb_usage.lstrip()
            return
        except SystemExit:
            return

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

        filename = args['<file>']
        system.save(filename)

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
            system.do_plot_obs(save=args['--save'])
        # plot residuals from fit
        if args['resid']:
            system.do_plot_resid(save=args['--save'])

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
                    system.do_plot_extras(e, save=args['--save'])
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
        
        verb = True if args['--verbose'] else False
        hf = float(args.pop('--hifac'))
        of = float(args.pop('--ofac'))
        fap = args['--fap']
        bfap = args['--bfap']
        showplot = not args['--noplot']

        # which periodogram should be calculated?
        per_fcn = None
        if args['--hoef']:
            per_fcn = periodograms.hoeffding
            name = 'Hoeffding'
        if args['--bgls']:
            per_fcn = periodograms.bgls
            name = 'Bayesian Generalized Lomb-Scargle'
        if args['--bayes']: 
            per_fcn = periodograms.bls
            name = 'Bayesian Lomb-Scargle'
        if args['--ls']: 
            per_fcn = periodograms.ls_PressRybicki
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
                if showplot:
                    system.per._plot(doFAP=fap, dobFAP=bfap, save=args['--save'])
            except AttributeError:
                system.per = per_fcn(system, hifac=hf, ofac=of)
                # system.per._output(verbose=verb)  # not ready
                if showplot:
                    system.per._plot(doFAP=fap, dobFAP=bfap, save=args['--save'])

        if args['bis']: # periodogram of the CCF's Bisector Inverse Slope
            system.bis_per = per_fcn(system, hifac=hf, ofac=of, quantity='bis')
            if showplot:
                system.bis_per._plot(doFAP=fap, dobFAP=bfap, save=args['--save'])

        if args['fwhm']: # periodogram of the CCF's fwhm
            system.fwhm_per = per_fcn(system, hifac=hf, ofac=of, quantity='fwhm')
            if showplot:
                system.fwhm_per._plot(doFAP=fap, dobFAP=bfap, save=args['--save'])

        if args['rhk']: # periodogram of rhk
            system.rhk_per = per_fcn(system, hifac=hf, ofac=of, quantity='rhk')
            if showplot:
                system.rhk_per._plot(doFAP=fap, dobFAP=bfap, save=args['--save'])

        if args['contrast']: # periodogram of contrast
            system.contrast_per = per_fcn(system, hifac=hf, ofac=of, quantity='contrast')
            if showplot:
                system.contrast_per._plot(doFAP=fap, dobFAP=bfap, save=args['--save'])

        if args['resid']: # periodogram of the residuals of the current fit
            system.resid_per = per_fcn(system, hifac=hf, ofac=of, quantity='resid')
            if showplot:
                system.resid_per._plot(doFAP=fap, dobFAP=bfap)

    @needs_local_scope
    @line_magic
    def clean(self, parameter_s='', local_ns=None):
        """ 
        Deconvolves the LS periodogram from the window function 
        using the CLEAN algorithm (Roberts et al. 1985)
        """
        # use default system or user defined
        try:
            if local_ns.has_key('default'):
                system = local_ns['default']
        except KeyError:
            from shell_colors import red
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        core.do_clean(system)

    @needs_local_scope
    @line_magic
    def wf(self, parameter_s='', local_ns=None):
        """ Calculate the spectral window function of the observations. 
        Type 'wf -h' for more help. """

        args = parse_arg_string('wf', parameter_s)
        if args == 1: return
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
            if system.model is None: system.model = {}
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
    def set_fit(self, parameter_s='', local_ns=None):
        from shell_colors import red
        # args = parse_arg_string('fit', parameter_s)
        # if args == 1: return
        #print args

        # verb = True if args['--verbose'] else False

        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        core.do_set_fit(system)


    @needs_local_scope
    @line_magic
    def correlate(self, parameter_s='', local_ns=None):
        from shell_colors import red
        args = parse_arg_string('correlate', parameter_s)
        if args == 1: return
        #print args

        verb = args['--verbose']
        rem = args['--remove']

        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        var1 = args['<var1>']
        var2 = args['<var2>']
        result = core.do_correlate(system, vars=(var1, var2), verbose=verb, remove=rem)

    @needs_local_scope
    @line_magic
    def de(self, parameter_s='', local_ns=None):
        """ Run the differential evolution algorithm minimization - stub """
        from shell_colors import red
        ## take care of arguments
        try:
            args = parse_arg_string('de', parameter_s)
        except DocoptExit:
            print de_usage.lstrip()
            return
        except SystemExit:
            return
        
        ngen = int(args.pop('--ngen'))
        npop = int(args.pop('--npop'))
        
        # default system?
        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        core.do_diffevol(system, npop=npop, ngen=ngen)
        system.do_plot_drift()
        system.do_plot_fit()

    @needs_local_scope
    @line_magic
    def demc(self, parameter_s='', local_ns=None):
        """ Run the Differential Evolution MCMC. - stub"""
        from shell_colors import red
        ## take care of arguments
        try:
            args = parse_arg_string('demc', parameter_s)
        except DocoptExit:
            print demc_usage.lstrip()
            return
        except SystemExit:
            return

        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return
        print args

        zfile = args.pop('<zipfile>')

        results = core.do_demc(system, zfile=zfile, burnin=0)
        return results
        # system.do_plot_fit()

    @needs_local_scope
    @line_magic
    def gen(self, parameter_s='', local_ns=None):
        """ Run the genetic algorithm minimization - stub """
        from shell_colors import red

        ## take care of arguments
        try:
            args = parse_arg_string('gen', parameter_s)
        except DocoptExit:
            print de_usage.lstrip()
            return
        except SystemExit:
            return
        
        ngen = int(args.pop('--ngen'))
        npop = int(args.pop('--npop'))


        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        core.do_genetic(system, npop=npop, ngen=ngen)
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
    def remove_rotation(self, parameter_s='', local_ns=None):
        """ Remove rotation period and harmonics """
        from shell_colors import red
        
        try:
            args = parse_arg_string('rrot', parameter_s)
        except DocoptExit:
            print rrot_usage.lstrip()
            return
        except SystemExit:
            return

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

        print args
        prot = args['--prot']
        nrem = int(args['--nrem'])
        fwhm = args['fwhm']
        rhk = args['rhk']

        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return
            
        core.do_remove_rotation(system, prot=prot, nrem=nrem, fwhm=fwhm, rhk=rhk)
        # core.do_genetic(system)
        # system.do_plot_fit()

    @needs_local_scope
    @line_magic
    def add_noise(self, parameter_s='', local_ns=None):
        try:
            args = parse_arg_string('add_noise', parameter_s)
        except DocoptExit:
            print addnoise_usage.lstrip()
            return
        except SystemExit:
            return

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

        # print args
        noise = float(args['<number>'])
        kms = system.error.mean() < 0.01
        if kms:
            system.error = sqrt(system.error**2 + (noise / 1000)**2)
        else:
            system.error = sqrt(system.error**2 + noise**2)

    @needs_local_scope
    @line_magic
    def lowpass(self, parameter_s='', local_ns=None):
        from shell_colors import blue
        # try:
        #     args = parse_arg_string('add_noise', parameter_s)
        # except DocoptExit:
        #     print addnoise_usage.lstrip()
        #     return
        # except SystemExit:
        #     return

        # use default system or user defined
        try:
            system = local_ns['default']
        except KeyError:
            from shell_colors import red
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return 

        core.do_lowpass_filter(system) 
        # f = 1./60

        # temp = lopast(system.extras.rhk, system.time, f)

        # print plt

        # system.vrad = system.vrad - temp


    @needs_local_scope
    @line_magic
    def killall(slef, parameter_s='', local_ns=None):
        from matplotlib.pyplot import close
        close('all')


    @needs_local_scope
    @line_magic
    def nest(self, parameter_s='', local_ns=None):
        """ Start the MultiNest analysis and handle data interaction and IO """
        from shell_colors import red

        try:
            args = parse_arg_string('nest', parameter_s)
        except DocoptExit:
            print nest_usage.lstrip()
            return
        except SystemExit:
            return
        # print args

        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        user = args['-u']
        resume = args['-r']
        verbose = args['-v']
        gp = args['--gp']
        jitter = args['--jitter']
        if gp and jitter:
            msg = red('ERROR: ') + '--gp and --jitter are incompatible'
            clogger.fatal(msg)
            return
        doplot = not args['--noplot']
        saveplot = args['--saveplot']
        dofeedback = args['--feed']
        doMAPfeedback = args['--MAPfeed']
        maxp = int(args['--maxp'])
        restart = args['--restart']
        nml_path = args['--nml']

        startp = args['--startp']
        if startp is not None:
            startp = [int(i) for i in startp.split(',')]
        else:
            startp = []

        try: 
            ncpu = int(args['--ncpu'])
        except TypeError:
            ncpu = None

        train_quantity = args['--train'] if bool(args['--train']) else None
        skip_train_mcmc = args['--skip-mcmc']
        lin_quantity = args['--lin'] if bool(args['--lin']) else None

        if bool(args['--train']) and not system.is_in_extras(train_quantity):
            msg = red('ERROR: ') + 'The name "%s" is not available in extras.\n' % train_quantity
            clogger.fatal(msg)
            return

        core.do_multinest(system, user, gp, jitter, maxp=maxp,
                          resume=resume, ncpu=ncpu, verbose=verbose,
                          training=train_quantity, skip_train_mcmc=skip_train_mcmc, lin=lin_quantity, 
                          doplot=doplot, saveplot=saveplot, feed=dofeedback, MAPfeed=doMAPfeedback,
                          restart=restart, nml=nml_path, startp=startp)

    @needs_local_scope
    @line_magic
    def dnest(self, parameter_s='', local_ns=None):
        pass

    @needs_local_scope
    @line_magic
    def dnest(self, parameter_s='', local_ns=None):
        from shell_colors import red

        if local_ns.has_key('default'):
            system = local_ns['default']
        else:
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        core.do_RJ_DNest3(system)

    @needs_local_scope
    @line_magic
    def restrict(self, parameter_s='', local_ns=None):
        """ Select data based on date, SNR or radial velocity accuracy.
        Type 'restrict -h' for more help
        """
        from shell_colors import yellow, blue, red
        args = parse_arg_string('restrict', parameter_s)

        if args == DocoptExit:
            msg = yellow('Warning: ') + "I'm not doing anything. Type restrict -h for help"
            clogger.fatal(msg)
            return

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
            core.do_restrict(system, 'gui')

        if args['--index']:
            core.do_restrict(system, 'index', args['--index'], noask=args['--noask'])

        # if args['--gui'] or args['--index']:
        #     if args['--index']:
        #         ind_to_remove = map(int, args['--index'].split(','))
        #         ind_to_remove = [i-1 for i in ind_to_remove]
        #         for i in ind_to_remove:
        #             x, y = take(system.time, i), take(system.vrad, i)
        #             msg = blue('INFO: ') + 'going to remove observation %d -> %8.2f, %8.2f\n' % (i+1, x, y)
        #             clogger.info(msg)
        #     else:
        #         ind_to_remove = selectable_plot(system, style='ro')

        #     n = len(ind_to_remove)
        #     if n == 0:
        #         msg = blue('    : ') + 'Not removing any observations'
        #         clogger.info(msg)
        #         return

        #     if args['--noask'] or ask_yes_no(red('    : ') + 'Are you sure you want to remove %d observations? (Y/n) ' % n, default=True):
        #         system.provenance.values()[0][1] += n
        #         # remove observations with indices ind_to_remove from
        #         # system.(time,vrad,error); leave *_full arrays intact
        #         system.time = delete(system.time, ind_to_remove)
        #         system.vrad = delete(system.vrad, ind_to_remove)
        #         system.error = delete(system.error, ind_to_remove)
        #         # remove observations with indices ind_to_remove from
        #         # system.extras.*; leave system.extras_full.* arrays intact
        #         for i, arr in enumerate(system.extras):
        #             field_name = system.extras._fields[i]
        #             replacer = {field_name:delete(arr, ind_to_remove)}
        #             system.extras = system.extras._replace(**replacer)
        #         msg = blue('    : ') + 'Done'
        #         clogger.info(msg)

        #         # delete system.per to force re-calculation
        #         try:
        #             del system.per
        #         except AttributeError:
        #             pass
        #     else:
        #         msg = blue('    : ') + 'Not removing any observations.'
        #         clogger.info(msg)                


    @needs_local_scope
    @line_magic
    def rotation(self, parameter_s='', local_ns=None):
        """ Calculate rotation period from activity-rotation relation"""
        args = parse_arg_string('rotation', parameter_s)
        # print args

        # use default system or user defined
        try:
            if local_ns.has_key('default') and not args['-n']:
                system = local_ns['default']
            else:
                system_name = args['SYSTEM']
                system = local_ns[system_name]
        except KeyError:
            from shell_colors import red
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        core.get_rotation_period(system)

    @needs_local_scope
    @line_magic
    def create(self, parameter_s='', local_ns=None):
        if '-h' in parameter_s:
            print create_usage
            return
        # print parameter_s
        
        if '--gui' in parameter_s:
            core.load_plugin('create_GUI')
            return

        core.do_create_planets(parameter_s)

    @needs_local_scope
    @line_magic
    def to_mps(self, parameter_s='', local_ns=None):
        # Convert data to meters per second, if in km/s
        from shell_colors import blue

        args = parse_arg_string('to_mps', parameter_s)
        # print args

        # use default system or user defined
        try:
            if local_ns.has_key('default') and not args['-n']:
                system = local_ns['default']
            else:
                system_name = args['SYSTEM']
                system = local_ns[system_name]
        except KeyError:
            from shell_colors import red
            msg = red('ERROR: ') + 'Set a default system or provide a system '+\
                                   'name with the -n option'
            clogger.fatal(msg)
            return

        if (min(system.error) < 0.01):
            msg = blue('INFO: ') + 'Converting to m/s'
            clogger.info(msg)

            system.vrad = (system.vrad - mean(system.vrad)) * 1e3
            system.error *= 1e3
            system.units = 'm/s'



    @needs_local_scope
    @line_magic
    def metalpoor(self, parameter_s='', local_ns=None):
        core.load_plugin('metalpoor')



def parse_arg_string(command, arg_string):
    """ Parse arguments for each of the commands. """
    # docopt does the heavy-lifting parsing, we just split the argument string
    # and catch the exceptions raised by -h or --help

    splitted = str(arg_string).split()

    if command is 'read':
        args = docopt(read_usage, splitted)

    if command is 'saverdb':
        args = docopt(saverdb_usage, splitted)        

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
            return DocoptExit

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

    if command is 'rrot':
        args = docopt(rrot_usage, splitted)

    if command is 'de':
        args = docopt(de_usage, splitted)

    if command is 'gen':
        args = docopt(gen_usage, splitted)

    if command is 'demc':
        args = docopt(demc_usage, splitted)        

    if command is 'add_noise':
        args = docopt(addnoise_usage, splitted)

    if command is 'rotation':
        args = docopt(rotation_usage, splitted)

    if command is 'to_mps':
        args = docopt(tp_mps_usage, splitted)

    if command is 'nest':
        args = docopt(nest_usage, splitted)

    return args
