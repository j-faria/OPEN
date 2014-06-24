#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
"""OPEN - software tools for exoplanet detection.

Usage:
  open.py
  open.py --nobanner
  open.py (-h | --help)
  open.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  --nobanner    Start OPEN without welcome message.
"""

# standard library imports
import sys
# intra-package imports
from OPEN.ipshell import EmbedShell, cfg
from OPEN.docopt import docopt
from OPEN._OPENversion import __version__

import numpy as np
import scipy as sp

if __name__ == "__main__":
    # version info
    version = 'OPEN v' + __version__
    # command-line argument parser
    args = docopt(__doc__, version=version)

    if '--art' in sys.argv:
        with open('./OPEN/asciiart.txt', 'r') as f: print f.read()

    banner =  u'Welcome to ' + version + '\n'
    banner += u'Created by João Faria | joao.faria@astro.up.pt\n\n'
    banner += 'Type "listcommands" for a list of the available commands\n'
    banner += 'For help on a specific command, type "command -h" or see http://j-faria.github.io/OPEN'
    ipyshell = EmbedShell(config=cfg, user_global_ns=globals(), user_ns={'np':np, 'sp':sp},
                          banner1=banner if (not args['--nobanner']) else '',
                          #banner2 = 'This is in banner 2',
                          exit_msg='Goodbye.')
    ipyshell.enable_pylab(import_all=False)
    ipyshell()