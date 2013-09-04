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


def raiseError(*args):

    if len(args) == 1:
        text = args[0]
        logger = None
    else:
        text = args[0]
        logger = args[1]

    if logger is not None:
        logger.write(cm.Fore.RED + '[ERROR]: ' + cm.Style.RESET_ALL +
                     text, newLine=True)
    else:
        print
        print cm.Fore.RED + '[ERROR]: ' + cm.Style.RESET_ALL + text

    sys.exit()

    return


def raiseWarning(*args, **kargs):

    if len(args) == 1:
        text = args[0]
        logger = None
    else:
        text = args[0]
        logger = args[1]

    doLog = kargs['doLog'] if 'doLog' in kargs else True
    doPrint = kargs['doPrint'] if 'doPrint' in kargs else True
    newLine = kargs['newLine'] if 'newLine' in kargs else True

    if logger is not None:
        logger.write(cm.Fore.YELLOW + '[Warning]: ' + cm.Style.RESET_ALL +
                     text, newLine=newLine, doLog=doLog, doPrint=doPrint)
    else:
        print cm.Fore.YELLOW + '[Warning]: ' + cm.Style.RESET_ALL + text

    return
