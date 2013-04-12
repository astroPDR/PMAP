#!/usr/bin/env python
# encoding: utf-8
"""

getRegions.py

Created by José Ramón Sánchez-Gallego on 2011-01-17.
Modified by José Ramón Sánchez-Gallego on 2011-07-04.
Copyright (c) 2011. All rights reserved.

This script is the main pipeline which reads the image names and options file and
performs all the necessary actions by calling the corresponding functions.

It works as follows:
(A) The user calls this pipeline with the FUV and HI images as arguments. An options file
    can be passed to the program. Otherwise, it will use the standard options.

"""

#JSH 2011-10-12: changed np.max and np.min to np.nanmax and np.nanmin respectively

import sys
import os
import re
import numpy as np
import pyfits as pf
import shutil as sh
from calcBackground import calcBackground
from clfind2dIDL import clfind2dIDL
from checkFUV import checkFUV
# from runCLFind import runCLFind

# Default options dictionary.
defOpts = dict(distance=None,
               log=True,
               logFile='test.log',
               verbose=False,
               fuvImage=None,
               hiImage=None,
               outputdir='./',
               fuvMaskFile='Default',
               fuvInteractive=True,
               fuvSigma=None,
               fuvhLevel=None,
               fuvnLevels=None,
               checkFUV=True,
               minDistance=500,
               scale=None,
               joinRegions=True,
               fluxContrast=3.0,
               peakContrast=3.0,
               alwaysRemoveClose=True,
               fuvMaskFileRej='Default',
               hiMaskFile='Default',
               hiInteractive=True,
               hiSigma=None,
               hihLevel=None,
               hinLevels=None)

# Dictionary with the final configuration options (either the defaults or
# the ones read from the config file)
configOpts = {}

# Fancy printing
lastLineStdOut = None
lastLineLog = None
unitLog = None
def fancyPrint(text, newLine=False, log=False, noPrint=False, flush=False):

    global lastLineStdOut, lastLineLog, unitLog
    textStdOut = textLog = text

    if newLine == True and flush == False and lastLineStdOut != '' and text!='':
        textStdOut = '\n' + str(text)

    if newLine == True and flush == False and lastLineLog != '' and text!='':
        textLog = '\n' + str(text)

    newLineChar = '\n'
    if flush: newLineChar = ''

    if noPrint == False:
        sys.stdout.write(str(textStdOut) + newLineChar)
        sys.stdout.flush()
        lastLineStdOut = textStdOut

    if log==True and configOpts.has_key('log') and configOpts['log']==True:
        if unitLog == None: unitLog = open(configOpts['logFile'], 'w')
        unitLog.write(textLog + newLineChar)
        lastLineStdOut = textLog
        unitLog.flush()

    return


# Opens the configuration file and retrieved the values. Returns a dictionary
def readConfigFile(fileOpts):
    configDic = {}
    lines = open(fileOpts, 'r').read().splitlines()
    for line in lines:
        strippedLine = line.strip()
        if strippedLine == '' or strippedLine[0] == '#': continue
        opt = strippedLine.split()[0]
        value = ' '.join(strippedLine.split()[1:])
        if '#' in value:
            comment = value[value.find('#'):].strip()
            value = value[0:value.find('#')].strip()
        if value == '':
            value = None
        elif value.lower()  in ['y', 'true', 'yes']:
            value = True
        elif value.lower()  in ['n', 'no', 'false']:
            value = False
        elif value.lower() in ['none', '', '-']:
            value = None
        else:
            try:
                value = float(value)
            except:
                pass
        configDic[opt] = value
    return configDic


# Parses the command-line call and returns the FUV and HI image names and the config file
# if any. In that case, it reads the config file. Otherwise, loads the default options.
def parseCommandLine():

    global configOpts

    import optparse

    parser = optparse.OptionParser()

    usage = 'usage: %prog [options] [FUV_Image] [HI_Image]'
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-c', '--config', dest='configFile',
        help='configuration file', default=None)
    parser.add_option('-y', '--default', dest='default',
        action='store_true', help='does not prompt. Uses default answers',
        default=False)

    (options, args) = parser.parse_args()

# Checks if the config file exists and loads it. Otherwise loads the default dir
    if options.configFile != None:
        if not os.path.exists(options.configFile):
            fancyPrint('Config file not found. Using defaults.', newLine=True)
            configFile = None
            configOpts = defOpts.copy()
        else:
            configFile = options.configFile
            configOpts = readConfigFile(configFile)
            fancyPrint('Configuration file %s loaded' % configFile,
                newLine=True)
    else:
        fancyPrint('Loading defaults options', newLine=True)
        configOpts = defOpts.copy()

    if os.path.exists(configOpts['logFile']): os.remove(configOpts['logFile'])

# The command-line images take precedence over the ones in the config file (in case
# they exist)
    if len(args) >= 1 and args[0] != 'None':
        configOpts['fuvImage'] = args[0]
    if len(args) >= 2 and args[1] != 'None':
        configOpts['hiImage'] = args[1]

# Checks for potential errors in the input
    if configOpts['fuvImage'] == None or configOpts['fuvImage'] == None:
        fancyPrint('FUV image or HI image not defined\n', newLine=True)
        sys.exit()

    if not os.path.exists(configOpts['fuvImage']) or not os.path.exists(configOpts['fuvImage']):
        fancyPrint('Some of the files don\'t exist\n', newLine=True)
        sys.exit()

    configOpts['default'] = options.default

    return configFile


# As parseCommandLine, but in this case it parses the call when getRegions is called as
# a subroutine
def parseInputData(inputData):

    global configOpts

    # Copies the input dictionary
    extraData = inputData.copy()

    # If the configuration file is defined in the inputData dictionary, it
    # loads it. Otherwise, uses the default options
    if 'configFile' in inputData:
        configFile = inputData['configFile']
        del extraData['configFile']
        try:
            configOpts = readConfigFile(configFile)
        except:
            fancyPrint('Config file not found. Using defaults.', newLine=True)
            configFile = None
            configOpts = defOpts.copy()
    else:
        configFile = None
        configOpts = defOpts.copy()

    # Sets the default option to false
    configOpts['default'] = False

    # Loads the remaining options in inputData
    for key,value in extraData.items():
        configOpts[key] = value

    # Checks for potential errors in the input
    if configOpts['fuvImage'] == None or configOpts['fuvImage'] == None:
        fancyPrint('FUV image or HI image not defined\n', newLine=True)
        sys.exit()

    if not os.path.exists(configOpts['fuvImage']) or not os.path.exists(configOpts['fuvImage']):
        fancyPrint('Some of the files don\'t exist\n', newLine=True)
        sys.exit()

    return configFile


def getFloat(text, default):
    value = None
    while type(value) not in [float, int]:
        value = raw_input(text)
        try:
            if value == '': return default
            value = float(value)
        except:
            pass
    return value


def checkCLFindOptions(root):
    everyThingOK = True
    try:
        if type(configOpts[root + 'Background']) not in [int, float]:
            everyThingOK = False
        if type(configOpts[root + 'hLevel']) not in [int, float, type(None)]:
            everyThingOK = False
        if type(configOpts[root + 'nLevels']) not in [int, float, type(None)]:
            everyThingOK = False
        if type(configOpts[root + 'Sigma']) not in [int, float]:
            everyThingOK = False
        if type(configOpts[root + 'MinPixels']) not in [int, float]:
            everyThingOK = False
    except:
        everyThingOK = False
    return everyThingOK


def calcParameters(root, sigmaSep = 3.):

    global configOpts

    image = configOpts[root + 'Image']
    sigma = configOpts[root + 'Sigma']
    background = configOpts[root + 'Background']

    lLevel = background * sigma

    configOpts[root + 'lLevel'] = lLevel

    if configOpts[root + 'hLevel'] in [None, 'None', 'N', '', '-1', -1.]:
        hLevel = np.nanmax(pf.getdata(image)) - background
        configOpts[root + 'hLevel'] = hLevel
    else:
        hLevel = configOpts[root + 'hLevel']

    if configOpts[root + 'nLevels'] in [None, 'None', 'N', '', '-1', -1.]:
        nLevels = np.floor((hLevel - lLevel)/(sigmaSep * background))
        configOpts[root + 'nLevels'] = nLevels
        levels = np.arange(lLevel, hLevel, sigmaSep * background)
    else:
        nLevels = configOpts[root + 'nLevels']
        levels = np.arange(lLevel, hLevel, (hLevel-lLevel)/nLevels)

    text = 'Image: %s\nBackground: %.3e\nSigma: %.1f\nLower level: %.3e\n'
    text += 'Higher level: %.3e\nNumber of levels: %d'
    fancyPrint('Running CLFind with the following parameters',
        newLine=True, log=True)
    fancyPrint(text % (image, background, sigma, lLevel, hLevel, nLevels),
        log=True)

    return levels


def getYN(text, returnOpt=None):
    global configOpts
    if configOpts['default']: return returnOpt
    opt = 'x'
    validKeys = {'y': True, 'n': False}
    if returnOpt != None: validKeys[''] = returnOpt
    while opt.lower() not in validKeys.keys():
        opt = raw_input(text)
    return validKeys[opt]



def getRegions(inputData):

    fancyPrint('')

    #####################################################################
    #######  This section reads the input and the config file ###########

    if inputData == None:
        configFile = parseCommandLine()
    else:
        configFile = parseInputData(inputData)

    fancyPrint('FUV Image: %s' % configOpts['fuvImage'], newLine = True, log=True)
    fancyPrint('HI Image: %s' % configOpts['hiImage'], log=True)
    fancyPrint('Config File: %s' % str(configFile), log=True)


    #####################################################################
    ##########################  FUV processing ##########################

    if configOpts['fuvInteractive'] == False:
        areCLParamsOK = checkCLFindOptions('fuv')
        if areCLParamsOK == False:
            fancyPrint('Some FUV parameters in config file are wrong\n',
                newLine=True, log=True)
            sys.exit()
        fancyPrint('FUV parameters loaded', newLine=True, log=True)

    else:

        fancyPrint('FUV: Running in interactive mode ... ', newLine=True, log=True)

        data = pf.getdata(configOpts['fuvImage'])
        fancyPrint('Image statistics', newLine=True, log=True)
        fancyPrint('Max: %9.3e' % np.nanmax(data), log=True)
        fancyPrint('Min: %9.3e' % np.nanmin(data), log=True)
        fancyPrint('Mean: %9.3e' % np.mean(data), log=True)
        fancyPrint('Median: %9.3e' % np.median(data), log=True)

        configOpts['fuvBackground'] = \
            calcBackground(configOpts['fuvImage'], fancyPrint)
        fancyPrint('Background: %.5e' % configOpts['fuvBackground'],
            newLine=True, log=False)

        fancyPrint('')
        configOpts['fuvSigma'] = getFloat('Sigma [3]: ', 3)
        configOpts['fuvhLevel'] = getFloat('Highest level [-1 for automatic]: ', -1)
        configOpts['fuvnLevels'] = getFloat('Number of levels  [-1 for automatic]: ', -1)
        configOpts['fuvMinPixels'] = getFloat('Minimum number of pixels [20]: ', 20)
        fancyPrint('')

    fuvLevels = calcParameters('fuv')

    if (configOpts['fuvMaskFile'].lower() == 'default') or (configOpts['fuvMaskFile'] == None):
        configOpts['fuvMaskFile'] = os.path.splitext(configOpts['fuvImage'])[0] + '_Msk.fits'

    doIDL = True
    if os.path.exists(configOpts['fuvMaskFile']):
        doIDL = getYN('\nFound {0}. Rerun CLFind? [y/N] '.format(configOpts['fuvMaskFile']),
            returnOpt=False)

    if doIDL:

        status, maskFile, idlLogFile = clfind2dIDL(configOpts['fuvImage'], fuvLevels,
            log=configOpts['log'], nPixMin=configOpts['fuvMinPixels'],
            verbose=configOpts['verbose'])

        if not status:
            fancyPrint('Problem found running IDL', newLine=True, log=True)
            fancyPrint('Kill IDL, check {0} and try again\n'.format(idlLogFile), log=True)
            sys.exit()

        sh.move(maskFile, configOpts['fuvMaskFile'])
        logData = open(idlLogFile, 'r').read().splitlines()
        sh.move(idlLogFile, 'clfind2dFUV.log')

    else:
        logData = open('clfind2dFUV.log', 'r').read().splitlines()

    for line in logData:
        fancyPrint(line, log=True, noPrint=True)
        if 'clumps found (' in line:
            m = re.match(r'(\d+)(.+)(\(+)(\d+)', line.replace(' ',''))
            nClumpsFUV = int(m.groups(0)[0])
            nClumpsRejFUV = nClumps = int(m.groups(0)[-1])

    fancyPrint('FUV mask file: {0}'.format(configOpts['fuvMaskFile']), newLine=True, log=True)
    fancyPrint('Number of clumps: {0}'.format(nClumpsFUV), log=True)
    fancyPrint('Number of clumps rejected: {0}'.format(nClumpsRejFUV), log=True)


    #####################################################################
    ##########################  FUV rejection ###########################

    if configOpts['checkFUV']:
        if (configOpts['fuvMaskFileRej'] == None) or (configOpts['fuvMaskFileRej'].lower() == 'default'):
            configOpts['fuvMaskFileRej'] = os.path.splitext(configOpts['fuvImage'])[0] + '_RejMsk.fits'

        doRejection = True
        if os.path.exists(configOpts['fuvMaskFileRej']) and doIDL == False:
            doRejection = getYN('\nFound {0}. Redo FUV rejection? [y/N] '.format(configOpts['fuvMaskFileRej']), returnOpt=False)

        if doRejection == True:
            try:
                checkFUV(configOpts, fancyPrint)
            except Exception as detail:
                raise detail
        else:
            fancyPrint('Not doing FUV rejection. Using image {0}'.format(configOpts['fuvMaskFileRej']),
                log=True, noPrint=True, newLine=True)


    # This is to avoid that ClumpFind runs on the HI image
    fancyPrint('')

    return configOpts

    #####################################################################
    ##########################  HI processing ###########################

    if configOpts['hiInteractive'] == False:
        areCLParamsOK = checkCLFindOptions('hi')
        if areCLParamsOK == False:
            fancyPrint('Some HI parameters in config file are wrong\n',
                newLine=True, log=True)
            sys.exit()
        fancyPrint('HI parameters loaded', newLine=True, log=True)

    else:

        fancyPrint('HI: Running in interactive mode ... ', newLine=True, log=True)

        data = pf.getdata(configOpts['hiImage'])
        fancyPrint('Image statistics', newLine=True, log=True)
        fancyPrint('Max: %9.3e' % np.nanmax(data), log=True)
        fancyPrint('Min: %9.3e' % np.nanmin(data), log=True)
        fancyPrint('Mean: %9.3e' % np.mean(data), log=True)
        fancyPrint('Median: %9.3e' % np.median(data), log=True)

        configOpts['hiBackground'] = \
            calcBackground(configOpts['hiImage'], fancyPrint)
        fancyPrint('Background: %.5e' % configOpts['hiBackground'],
            newLine=True, log=False)

        fancyPrint('')
        configOpts['hiSigma'] = getFloat('Sigma [3]: ', 3)
        configOpts['hihLevel'] = getFloat('Highest level [-1 for automatic]: ', -1)
        configOpts['hinLevels'] = getFloat('Number of levels  [-1 for automatic]: ', -1)
        configOpts['hiMinPixels'] = getFloat('Minimum number of pixels [20]: ', 20)
        fancyPrint('')

    hiLevels = calcParameters('hi')

    if (configOpts['hiMaskFile'].lower() == 'default') or (configOpts['hiMaskFile'] == None):
        configOpts['hiMaskFile'] = os.path.splitext(configOpts['hiImage'])[0] + '_Msk.fits'

    doIDL = True
    if os.path.exists(configOpts['hiMaskFile']):
        doIDL = getYN('\nFound {0}. Rerun CLFind? [y/N] '.format(configOpts['hiMaskFile']),
            returnOpt=False)

    if doIDL:

        status, maskFile, idlLogFile = clfind2dIDL(configOpts['hiImage'], hiLevels,
            log=configOpts['log'], nPixMin=configOpts['hiMinPixels'],
            verbose=configOpts['verbose'])

        if not status:
            fancyPrint('Problem found running IDL', newLine=True, log=True)
            fancyPrint('Kill IDL, check {0} and try again\n'.format(idlLogFile), log=True)
            sys.exit()

        sh.move(maskFile, configOpts['hiMaskFile'])
        logData = open(idlLogFile, 'r').read().splitlines()
        sh.move(idlLogFile, 'clfind2dHI.log')

    else:
        logData = open('clfind2dHI.log', 'r').read().splitlines()

    for line in logData:
        fancyPrint(line, log=True, noPrint=True)
        if 'clumps found (' in line:
            m = re.match(r'(\d+)(.+)(\(+)(\d+)', line.replace(' ',''))
            nClumpsHI = int(m.groups(0)[0])
            nClumpsRejHI = nClumps = int(m.groups(0)[-1])

    fancyPrint('HI mask file: {0}'.format(configOpts['hiMaskFile']), newLine=True, log=True)
    fancyPrint('Number of clumps: {0}'.format(nClumpsHI), log=True)
    fancyPrint('Number of clumps rejected: {0}'.format(nClumpsRejHI), log=True)


    # #####################################################################
    # ##########################  HI rejection ###########################
    #
    # if configOpts['checkHI']:
    #     if (configOpts['fuvMaskFileRej'] == None) or (configOpts['fuvMaskFileRej'].lower() == 'default'):
    #         configOpts['fuvMaskFileRej'] = os.path.splitext(configOpts['fuvImage'])[0] + '_RejMsk.fits'
    #
    #     doRejection = True
    #     if os.path.exists(configOpts['fuvMaskFileRej']) and doIDL == False:
    #         doRejection = getYN('\nFound {0}. Redo FUV rejection? [y/N] '.format(configOpts['fuvMaskFileRej']), returnOpt=False)
    #
    #     if doRejection == True:
    #         try:
    #             checkFUV(configOpts, fancyPrint)
    #         except Exception as detail:
    #             raise detail
    #     else:
    #         fancyPrint('Not doing FUV rejection. Using image {0}'.format(configOpts['fuvMaskFileRej']),
    #             log=True, noPrint=True, newLine=True)
    #
    #
    # #####################################################################

    fancyPrint('')

    return configOpts


if __name__ == '__main__':

    getRegions(None)

