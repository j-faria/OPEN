# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
"""
Provides support for printing colors in the terminal by using ANSI escape
sequences.

To use this code, you can do something like:
    print scolor.WARNING + "This is printed in yellow" + scolor.ENDC
or
    import shell_colors as sc
    print sc.red('This is printed in red')
    
In principle this will work on unix, linux (including macOS), and windows 
(provided you enable ansi.sys)
"""

class scolor:
    HEADER = '\033[95m'     # magenta
    OKBLUE = '\033[94m'     # blue
    OKGREEN = '\033[92m'    # green
    WARNING = '\033[93m'    # yellow
    FAIL = '\033[91m'       # red
    ENDC = '\033[0m'        # termination sequence 

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        
def magenta(message):
    """ Return a magenta(ed) message """
    return scolor.HEADER + message + scolor.ENDC

def blue(message):
    """ Return a blue(ed) message """
    return scolor.OKBLUE + message + scolor.ENDC
    
def green(message):
    """ Return a green(ed) message """
    return scolor.OKGREEN + message + scolor.ENDC
    
def yellow(message):
    """ Return a yellow(ed) message """
    return scolor.WARNING + message + scolor.ENDC
    
def red(message):
    """ Return a red(ed) message """
    return scolor.FAIL + message + scolor.ENDC
    
def clean(string):
    """ Clean a colored skin by removing escape sequences """
    return string[5:-4]
    
    
    
    
