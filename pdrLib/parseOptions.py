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
import numpy as np
from Error import raiseWarning, raiseError


class ConfigOptions(object):

    def __init__(self, configFile, createConfig=False, interactiveConfig=False,
                 verbose=True):

        self.configFile = configFile
        self.verbose = verbose
        self.defaultFile = os.path.dirname(__file__) + '/defaultOptions.dat'
        self.mandatoryOptions = os.path.dirname(__file__) + '/mandatoryOptions.dat'

        if createConfig is True:
            if interactiveConfig is False:
                # If we are creating a default configuration file
                self.creationStatus = self.createConfigFile()
            else:
                # If we are creating a configuration file in which the user interactively
                # inputs the values to use
                self.creationStatus = self.createConfigFileInteractive()
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
                    _config.set(_section, 'root',
                                os.path.splitext(os.path.basename(self.configFile))[0])

            for _item, _value in _config.items(_section):
                try:
                    _optionsDict[_item] = eval(_value)  # If possible, evaluates the data
                except:
                    _optionsDict[_item] = _value

                # If the field is empty (or starts with a comment), changes the value to None
                if isinstance(_optionsDict[_item], str):
                    if _optionsDict[_item].strip() == '' or \
                            _optionsDict[_item][0] in ['', ';', '#']:
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
            self.OptionsError('Some mandatory parameters are missing. Exiting now.')

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

    def clearComments(self, config):
        """
        Clears all the comments from a ConfigParser instance
        """

        _sections = config.sections()

        for _section in _sections:
            _fields = config.options(_section)
            for _field in _fields:
                _tmpValue = config.get(_section, _field).strip()
                if _tmpValue[0] == ';' or _tmpValue[0] == '#':
                    _tmpValue = ''
                config.set(_section, _field, _tmpValue)

        return

    def requestInput(self, text, default=None):
        """
        Recursively request an input until an acceptable
        value is introduced. If default is not None and
        and empty value is input, it returns default.
        """

        while True:
            if default is None:
                _defaultStr = ''
            else:
                _defaultStr = ' [%s]' % str(default)

            try:
                _answer = raw_input('%s%s: ' % (text, _defaultStr))

                if _answer.strip() == '':
                    if default is None:
                        pass
                    else:
                        return default
                else:
                    return _answer.strip()
            except KeyboardInterrupt:
                sys.exit(0)
            finally:
                pass

    def getGalParams(self, galaxy):
        """
        Uses astroquery and Simbad to get the galaxy parameters.
        If fails, returns 0 if astroquery is not presents, 1 if the version
        is too old or 2 if something went wrong when parsing Simbad.
        Otherwise, returns a tuple (distance, pa, incl, c_RA, c_DEC, R25)
        """

        from distutils.version import LooseVersion

        import warnings
        warnings.simplefilter('ignore')
        try:
            # This is only assured to work with astroquery >= 0.2
            # so we check to make sure that it is installed
            import astroquery
            _versionAQ = astroquery.__version__
            if LooseVersion(_versionAQ) < LooseVersion('0.2'):
                return 1
        except:
            # If there is no astroquery installed at all
            return 0

        from astroquery.simbad import Simbad

        Simbad.add_votable_fields('dim', 'rvel')
        try:
            _table = Simbad.query_object(galaxy)
        except:
            return 2

        # If no data is found for the galaxy
        if len(_table) == 0:
            return 2

        _ra = _table['RA'][0]
        try:
            _raSplit = _ra.split()
            _raHMS = '%sh%sm%ss' % (_raSplit[0], _raSplit[1], _raSplit[2])
        except:
            _raHMS = None

        _dec = _table['DEC'][0]
        try:
            _decSplit = _dec.split()
            _decDMS = '%sd%sm%ss' % (_decSplit[0], _decSplit[1], _decSplit[2])
        except:
            _decDMS = None

        _pa = _table['GALDIM_ANGLE'][0]
        if _pa == '':
            _pa = None

        _incl = _table['GALDIM_INCL'][0]
        if _incl == '':
            _incl = None

        _incl = _table['GALDIM_INCL'][0]
        if _incl == '':
            _incl = None

        _vrec = _table['RVel_Rvel'][0]
        if _vrec == '':
            _vrec = None
            _distance = None
            _r25Kpc = None
        else:
            _distance = _vrec / 71.   # Hubble's law

            # Calculates R25 in kpc from the major axis of the galaxy
            _majAxis = _table['GALDIM_MAJAXIS'][0]
            if _majAxis == '':
                _r25Kpc = None
            else:
                _r25Kpc = _distance * np.arctan(np.deg2rad(_majAxis / 2. / 60.)) * 1e3

        return (_distance, _pa, _incl, _raHMS, _decDMS, _r25Kpc)

    def createConfigFileInteractive(self):
        """
        Creates a configuration file after asking the user to input the key
        parameter via stdin.
        """

        from calcScale import calcScale

        # First, we create an empty configuration file
        _resultCreateConfig = self.createConfigFile()
        if _resultCreateConfig is False:
            raiseError('Impossible to create configuration file')

        # Reads the configuration file
        _config = ConfigParser.SafeConfigParser()
        _config.optionxform = str
        _config.read(self.configFile)

        # self.clearComments(_config)

        print(cm.Fore.YELLOW + '\nThis process will ask you the minimum parameter needed to run')
        print('the pipeline, and will generate a configuration file using')
        print('those values.\n' + cm.Style.RESET_ALL)

        # Filenames
        _defaultRoot = os.path.splitext(os.path.basename(self.configFile))[0]
        _root = self.requestInput('root', default=_defaultRoot)
        _fileNamesFieldsToChange = ['data_uv', 'data_hi', 'logFile',
                                    'data_results_unfiltered', 'data_results']

        _config.set('Logging options', 'root', _root)

        for _fileNameFieldsToChange in _fileNamesFieldsToChange:
            _tmpValue = _config.get('Logging options', _fileNameFieldsToChange)
            _tmpValue = _tmpValue.replace(_root, '%(root)s')
            _config.set('Logging options', _fileNameFieldsToChange, _tmpValue)

        # Pixel scales
        _imageNameFields = [['fuvImage', 'FUV_Pixel_Scale'],
                           ['hiImage', 'HI_Pixel_Scale']]
        for _imageNameField, _scaleNameField in _imageNameFields:
            _image = self.requestInput(_imageNameField)
            _config.set('Images', _imageNameField, _image)
            try:
                _scaleFromImage = calcScale(_image)
            except:
                _scaleFromImage = None
            _scale = self.requestInput(_scaleNameField, default=_scaleFromImage)
            _config.set('Images', _scaleNameField, str(_scale))

        # Conversion factors
        print(cm.Fore.YELLOW + '\nInput the following conversion factors for the images.')
        print ('FUV: CPS to erg sec-1 cm-2 A-1')
        print ('HI: CPS to cm-2')
        print('The default values proposed are the conversion factor for GALEX and THINGS')
        print ('respectively. Use with caution.' + cm.Style.RESET_ALL)

        _fuvScale = self.requestInput('FUVscale', default='1.40e-15')
        _config.set('Images', 'FUVscale', _fuvScale)
        _hiScale = self.requestInput('HIscale', default='4.907e19')
        _config.set('Images', 'HIscale', _hiScale)

        # Galaxy parameters (tries to use astroquery and Simbad)
        print(cm.Fore.YELLOW + '\nI\'ll try to connect to Simbad and download the parameters')
        print('of the galaxy.' + cm.Style.RESET_ALL)

        # Asks for the name of the galaxy (in case is different form the filaname root)
        _galName = self.requestInput('Galaxy name', default=_root)

        # Tries to run astroquery. If something fails, _galParamsReturn=2 (unknown error)
        try:
            _galParamsReturn = self.getGalParams(_galName)
        except:
            _galParamsReturn = 2

        # Depending on the flag received from self.getGalParams
        sys.stdout.write(cm.Fore.YELLOW)
        sys.stdout.flush()
        if _galParamsReturn == 0:
            # No data obtained. Makes _galParams = None so that the user needs
            # to input all the values
            print ('Astroquery is not installed. No parameters downloaded.')
            _galParams = (None, None, None, None, None, None)
        elif _galParamsReturn == 1:
            print ('Astroquery is present but a version >= 0.2 is required.')
            _galParams = (None, None, None, None, None, None)
        elif _galParamsReturn == 2:
            print ('There was an unknown error and the parameters could not be downloaded.')
            _galParams = (None, None, None, None, None, None)
        else:
            # If _galParamsReturn is not and int, it means that getGalParams has returned
            #  a tuple with the parameters. We make _galParams = _galParamsReturn
            print('The parameters were successfully downloaded and will')
            print('be shown as default options.')
            _galParams = _galParamsReturn
        print(cm.Style.RESET_ALL)

        _galFields = ['distance', 'pa', 'incl', 'c_RA', 'c_DEC', 'R25']
        for _nn, _galField in enumerate(_galFields):
            _tmpValue = self.requestInput(_galField, default=_galParams[_nn])
            _config.set('Physical properties', _galField, str(_tmpValue))

        _tmpValue = self.requestInput('ext', default=None)
        _config.set('Physical properties', 'ext', str(_tmpValue))

        print('')

        # Background levels
        _tmpFUVSigma = _config.get('FUV', 'fuvSigma')
        if _tmpFUVSigma.lower() == 'none':
            _fuvSigma = self.requestInput('fuvSigma', default=None)
        else:
            _fuvSigma = self.requestInput('fuvSigma', default=_tmpFUVSigma)
        _config.set('FUV', 'fuvSigma', _fuvSigma)

        print('')
        # Dust-to-gas-ratio

        # Before removing the default dust-gas section (which contains the three method)
        # we get some values
        _dusterr = _config.get('Dust-to-gas ratio', 'dusterr')
        _slopeerror = _config.get('Dust-to-gas ratio', 'slopeerror')
        _solar = _config.get('Dust-to-gas ratio', 'solar')

        # Resets section
        _config.remove_section('Dust-to-gas ratio')
        _config.add_section('Dust-to-gas ratio')

        _typeDust = self.requestInput('dusttype (slope | map | const)', default='const')
        if _typeDust == 'const':
            _dustVal = self.requestInput('dustval')
            _config.set('Dust-to-gas ratio', 'dusttype', _typeDust)
            _config.set('Dust-to-gas ratio', 'dustval', str(_dustVal))
            _config.set('Dust-to-gas ratio', 'dusterr', _dusterr)
        elif _typeDust == 'slope':
            _dustSlope = self.requestInput('slope')
            _dustOffset = self.requestInput('offset')
            _config.set('Dust-to-gas ratio', 'dusttype', _typeDust)
            _config.set('Dust-to-gas ratio', 'slope', str(_dustSlope))
            _config.set('Dust-to-gas ratio', 'offset', str(_dustOffset))
            _config.set('Dust-to-gas ratio', 'slopeerror', _slopeerror)
        elif _typeDust == 'map':
            _dustMap = self.requestInput('mapname')
            _dustMapError = self.requestInput('maperror')
            _config.set('Dust-to-gas ratio', 'dusttype', _typeDust)
            _config.set('Dust-to-gas ratio', 'mapname', str(_dustMap))
            _config.set('Dust-to-gas ratio', 'maperror', str(_dustMapError))
            _config.set('Dust-to-gas ratio', 'solar', _solar)

        with open(self.configFile, 'wb') as _unitConfigFile:
            _config.write(_unitConfigFile)

        print(cm.Fore.YELLOW + '\nThe configuration file %s has been generated.' % self.configFile)
        print('Please, check the values before running the pipeline.')
        print('The complete set of parameters and their detailed explanation ')
        print('can be found in https://github.com/astroPDR/PMAP.' + cm.Style.RESET_ALL)

        return True

    def logOpts(self, logger):
        """
        Saves the configuration options to the logger
        """

        for _option in sorted(self.options.keys(), key=str.lower):
            logger.write('%s = %s' % (_option, str(self.options[_option])), doPrint=False)

        return

    def OptionsError(self, error, exit=True):
        print cm.Fore.RED + '[ERROR]: ' + cm.Style.RESET_ALL + error
        if exit is True:
            sys.exit(2)
