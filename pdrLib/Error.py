#!/usr/bin/env python
# encoding: utf-8
"""
Error.py

Created by José Sánchez-Gallego on 15 Apr 2013.
Copyright (c) 2013. All rights reserved.

A series of routines to manage exceptions

"""

import sys
import colorama as cm


def raiseError(text, logger):

    logger.write(cm.Fore.RED + '[ERROR]: ' + cm.Style.RESET_ALL +
                 text, newLine=True)
    sys.exit()

    return
