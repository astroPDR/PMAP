#!/usr/bin/env python
# encoding: utf-8
"""
parseOptions.py

Created by José Sánchez-Gallego on 13 Apr 2013.
Copyright (c) 2013. All rights reserved.

This class creates an optiond instance with the configuration
parameters.

"""

import os
import ConfigParser
import warnings
import colorama as cm
import sys
from Error import raiseWarning


class ConfigOptions(object):

    def __init__(self, configFile, createConfig=False, verbose=True):

        self.configFile = configFile
        self.verbose = verbose
        self.defaultFile = os.path.dirname(__file__) + '/defaultOptions.dat'
        self.mandatoryOptions = os.path.dirname(__file__) + '/mandatoryOptions.dat'

        if createConfig is True:
            self.creationStatus = self.createConfigFile()
        else:
            self.creationStatus = False
            self.defaultOpts = self.parseConfigFile(self.defaultFile, dust=False)
            self.options = self.parseConfigFile(self.configFile)
            self.checkOptions()

        return

    def parseConfigFile(self, fileToParse, dust=True):
        """
        Reads the configuration file and returns a dictionary with all the options.
        """

        if not os.path.exists(fileToParse):
            raise IOError('Config file %s not found' % fileToParse)

        _config = ConfigParser.SafeConfigParser()
        _config.optionxform = str
        _config.read(fileToParse)

        _sections = _config.sections()
        _optionsDict = {}
        for _section in _sections:

            if _section == 'Dust-to-gas ratio' and dust is False:
                continue

            if (fileToParse == self.configFile) and ('logging' in _section.strip().lower()):
                if _config.has_option(_section, 'root') and \
                        _config.get(_section, 'root').lower() == 'default':
                    _config.set(_section, 'root', os.path.splitext(os.path.basename(self.configFile))[0])

            for _item, _value in _config.items(_section):
                try:
                    _optionsDict[_item] = eval(_value)  # If possible, evaluates the data
                except:
                    _optionsDict[_item] = _value

                # If the field is empty (or starts with a comment), changes the value to None
                if type(_optionsDict[_item]) is str and _optionsDict[_item][0] in [';', '']:
                    _optionsDict[_item] = None
                # This is to avoid eval from converting map into a function
                if _value.lower() == 'map':
                    _optionsDict[_item] = _value

        return _optionsDict

    def checkOptions(self):
        """
        Makes sure that all the necessary options are present. If a non-essential option is
        not present, used the default value. Raises warnings or error depending on the missing
        value.
        """

        _exitAtTheEnd = False
        _errorOpts = ConfigParser.SafeConfigParser()
        _errorOpts.optionxform = str
        _errorOpts.read(self.mandatoryOptions)

        _error = _errorOpts.get('Error', 'errorOpts').replace(' ', '').split(',')
        _warning = _errorOpts.get('Warning', 'warningOpts').replace(' ', '').split(',')

        if 'dusttype' not in self.options:
            self.OptionsError('dusttype not defined', exit=False)
            _exitAtTheEnd = True

        else:
            _error.append(_errorOpts.get('Dust',
                                         self.options['dusttype']).replace(' ', '').split(','))

        _fields = self.defaultOpts.keys()
        for _field in sorted(_fields, key=str.lower):

            if _field not in self.options.keys() and self.defaultOpts[_field] is not None \
                    and self.verbose is True:
                self.options[_field] = self.defaultOpts[_field]
                warnings.warn('Option %s not present. Using default value \'%s\'.' % (_field,
                              str(self.defaultOpts[_field])))

            if self.options[_field] is None and _field in _error:
                self.OptionsError('Option %s is mandatory' % _field, exit=False)
                _exitAtTheEnd = True
            elif self.options[_field] is None and _field in _warning and self.verbose is True:
                # warnings.warn(cm.Fore.YELLOW + '[WARNING]: ' + cm.Style.RESET_ALL +
                raiseWarning('Option %s is not defined. ' % _field +
                             'This won\'t stop the pipeline but may or may not be a problem.')

        if _exitAtTheEnd is True:
            self.OptionsError('Some mandatory parameters are missing. Exitting now.')

        return

    def createConfigFile(self):
        """
        Creates a copy of the sample configuration file
        """

        import shutil

        try:
            shutil.copyfile(self.defaultFile, self.configFile)
            return True
        except:
            return False

    def logOpts(self, logger):
        """
        Saves the config options to the logger
        """

        for _option in sorted(self.options.keys(), key=str.lower):
            logger.write('%s = %s' % (_option, str(self.options[_option])), doPrint=False)

        return

    def OptionsError(self, error, exit=True):
        print cm.Fore.RED + '[ERROR]: ' + cm.Style.RESET_ALL + error
        if exit is True:
            sys.exit(2)
