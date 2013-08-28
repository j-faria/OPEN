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
  --speed=<kn>  Speed in knots [default: 10].
  --moored      Moored (anchored) mine.
  --drifting    Drifting mine.
  --nobanner    Start OPEN without welcome message.
"""

import sys
from ipshell import EmbedShell, cfg
from docopt import docopt

from default import default_a

if __name__ == "__main__":

    # version info
    from _version import __version__
    version = 'OPEN v' + __version__    
    # command-line argument parser
    args = docopt(__doc__, version=version)
    print args
    
    if '--art' in sys.argv:
        with open('./OPEN/asciiart.txt', 'r') as f: print f.read()
    
    banner =  u'Welcome to OPEN v0.0 \n'
    banner += u'Copyright (c) João Faria | joao.faria@astro.up.pt\n\n'
    banner +=  'Type "listcommands" for a list of the available commands\n'
    banner +=  'For help on a specific command, type "command?" or see www.'
    ipyshell = EmbedShell(config=cfg, user_global_ns=globals(),
                          banner1 = banner if (not args['--nobanner']) else '',
                          #banner2 = 'This is in banner 2',
                          exit_msg = 'Goodbye.')
    ipyshell()

