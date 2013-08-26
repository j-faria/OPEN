# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#
import logging


# The following code extends Python's logging module with a custom log level
# called VERBOSE. 
# This log level sits between the predefined INFO and DEBUG levels and is used 
# to convey "detailed information that should be understandable to experienced 
# users to provide insight in the software's behavior; a sort of high level 
# debugging information".
# According the logging documentation this is not ideal but I found OPEN did in
# fact need this aditional level.

# Define the numeric value and label of the verbose level.
logging.VERBOSE = 15
logging.addLevelName(logging.VERBOSE, 'VERBOSE')

class VerboseLogger(logging.Logger):
    """
    Custom logger class that supports the additional logging level "verbose" 
    whose severity sits between "info" and "debug".
    """
    def __init__(self, *args, **kw):
        """
        Initialize the superclass and define the custom "verbose" log level.
        """
        logging.Logger.__init__(self, *args, **kw)

    def verbose(self, *args, **kw):
        """
        Log a verbose message: A message that we would like to log in verbose
        mode (-v) as a sort of high level debugging information (whereas
        logger.debug() is used to log low level information). This method has
        the same contract as the existing methods for logging a message.
        """
        self.log(logging.VERBOSE, *args, **kw)
        

# instantiate our custom logger at the default level of INFO
# every module should use this logger
clogger = VerboseLogger(__name__)
clogger.addHandler(logging.StreamHandler())
clogger.setLevel(logging.INFO)

