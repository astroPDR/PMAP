#!/usr/bin/env python
# encoding: utf-8
"""
fancyPrint.py

Created by José Sánchez-Gallego on 14 Apr 2013.
Copyright (c) 2013. All rights reserved.

This class createas a fancyPrint instance for personalised logging.

"""

import os
import re
import sys


class fancyPrint(object):

    def __init__(self, writeLog=True, logFile=None, verbose=False):

        self.logText = '\n'
        self.printText = '\n'
        self.verbose = verbose
        self.logFile = logFile
        self.writeLog = writeLog

        if self.writeLog is True and self.logFile is not None:
            self.unitOut = open(self.logFile, 'w')

        return

    def write(self, text, newLine=False, doLog=True, doPrint=None):

        if doPrint is None:
            doPrint = self.verbose

        ansiEscape = re.compile(r'\x1b[^m]*m')

        if newLine is True:
            if self.printText == '\n' or self.printText[-2:] != '\n\n':
                textToPrint = '\n' + text
            else:
                textToPrint = text
            if self.logText >= '\n' or self.logText[-2:] != '\n\n':
                textToLog = '\n' + ansiEscape.sub('', text)
            else:
                textToLog = ansiEscape.sub('', text)
        else:
            textToPrint = text
            textToLog = ansiEscape.sub('', text)

        if doPrint is True:
            sys.stdout.write(textToPrint + '\n')
            sys.stdout.flush()
            self.printText += textToPrint + '\n'

        if self.writeLog is True and doLog is True:
            self.unitOut.write(textToLog + '\n')
            self.unitOut.flush()
            os.fsync(self.unitOut.fileno())
            self.logText += textToLog + '\n'

        return
