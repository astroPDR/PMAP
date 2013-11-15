#! /usr/bin/env python
# encoding: utf-8
"""
PMAP.py

This file is the main routine for the PMAP pipeline.
The pipeline can be called using the command

python PMAP.py <configFile>

A new, empty configuration file can be created using

python PMAP.py -c <configFile>

or a new configuration file can be created interactively
with

python PMAP.py -i <configFile>

In the latter case, the pipeline will prompt for the
minimum needed values to run the code.

Create by Jonathan Heiner and José Sánchez-Gallego.

"""

import pdrLib as pdr
import os
import numpy as np
from pdrLib.Error import raiseError


def savePixScales(configOpts):

    # Calculates the FUV pixel scale
    if configOpts['FUV_Pix_Scale'] is None:
        configOpts['FUV_Pix_Scale'] = pdr.calcScale(configOpts['fuvImage'])
    else:
        if not isinstance(configOpts['FUV_Pix_Scale'], (float, int)):
            raiseError('FUV pixel scale is not a float.')

    # FUV pixel scale in pc
    configOpts['FUV_PcPerPixel'] = 2e6 * configOpts['distance'] * \
        np.tan(np.deg2rad(configOpts['FUV_Pix_Scale'] / 3600. / 2.))

    # Idem for the HI one
    if configOpts['HI_Pix_Scale'] is None:
        configOpts['HI_Pix_Scale'] = pdr.calcScale(configOpts['hiImage'])
    else:
        if not isinstance(configOpts['HI_Pix_Scale'], (float, int)):
            raiseError('HI pixel scale is not a float.')

    # HI pixel scale in pc
    configOpts['HI_PcPerPixel'] = 2e6 * configOpts['distance'] * \
        np.tan(np.deg2rad(configOpts['HI_Pix_Scale'] / 3600. / 2.))

    return


def checkAstropyVersion(configOpts):
    """
    Checks the version of astropy installed (if any). Adds two
    options to configOpts:
        - apVersion is the real version of astropy from astropy.__version__
        - apVersionSimple is a simplified version that can be 0.2.4 or 0.3
            as the syntaxis frequently changes from one to the other.
    """

    from distutils.version import LooseVersion

    try:
        import astropy
        isAstropyInstalled = True
    except:
        isAstropyInstalled = False

    if isAstropyInstalled is False:
        raiseError('Astropy cannot be imported. Please install it and try again.')

    apVersion = astropy.__version__
    if LooseVersion(apVersion) < LooseVersion('0.2.4'):
        raiseError('Astropy version is < 0.2.4. PMAP cannot continue until you update.')
    elif (LooseVersion(apVersion) >= LooseVersion('0.2.4')) and \
            (LooseVersion(apVersion) < LooseVersion('0.3')):
        apVersionSimple = '0.2.4'
    else:
        apVersionSimple = '0.3'

    configOpts['apVersion'] = apVersion
    configOpts['apVersionSimple'] = apVersionSimple

    return


###############################################################################
###                                 MAIN                                    ###
###############################################################################

def PMAP(configFile, createConfig=False, interactiveConfig=False,
         verbose=True, overwrite=False):
    """
    This main routine sequentially calls all the subroutines needed for the PDR
    method to work.
    """

    # configOptionsInstance is a instance of the class ConfigOptions. For convenience,
    # configOpts allows a quick access to the configOptionsInstance.options dict.
    configOptionsInstance = pdr.ConfigOptions(configFile, createConfig=createConfig,
                                              interactiveConfig=interactiveConfig,
                                              verbose=verbose)

    if createConfig is True:
        if configOptionsInstance.creationStatus is True:
            return
        else:
            raise IOError('Impossible to create %s' % configFile)
    else:
        configOpts = configOptionsInstance.options

    # Creates the logger

    if configOpts['log'] is True:
        logger = pdr.fancyPrint(logFile=configOpts['logFile'], writeLog=True,
                                verbose=verbose)
    else:
        logger = pdr.fancyPrint(logFile=None, writeLog=False, verbose=verbose)

    logger.write('Pipeline run with options:\n', doPrint=False)
    configOptionsInstance.logOpts(logger)

    # Adds the verbose and overwrite options to the configOpts dictionary
    configOpts['verbose'] = verbose
    if 'overwrite' not in configOpts:
        configOpts['overwrite'] = overwrite
    else:
        if configOpts['overwrite'] is not True:
            configOpts['overwrite'] = overwrite
        else:
            pass

    # Calculates the scale of the HI and FUV images, for later use
    savePixScales(configOpts)

    # Checks astropy version
    checkAstropyVersion(configOpts)

    # Checks the FOV of the FUV and HI images and matchs them
    matchFOV = configOpts['matchFOV']
    if matchFOV is True:
        fov = pdr.matchFOV(configOpts['fuvImage'], configOpts['hiImage'])
    else:
        fov = None

    # Calls the getRegions routine, which produces the region masks for the FUV and HI images
    pdr.getRegions(configOpts, logger, fov=fov)

    #Todo: call reject regions on _Peaks.dat if desired

    fluxFUVTable = pdr.getUVFlux(configOpts, logger)
    hiData = pdr.getHIFlux(configOpts, fluxFUVTable, logger, fov=fov)

    # Now we apply the PDR mehod to the FUV and HI measurements
    pdrData = pdr.applyPDRMethod(configOpts, fluxFUVTable, hiData, logger)

    return pdrData


if __name__ == '__main__':
    """
    Reads command line arguments
    """

    from optparse import OptionParser

    usage = 'usage: python %prog configFile [runs pipeline]' + \
            '\n       python %prog -c [configFile] [generates a sample config file]' + \
            '\n       python %prog -i [configFile] [generates a sample config file interactively]'

    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--config', dest='sampleFile',
                      help='generates a default configuration file.',
                      default=False, action='store_true')
    parser.add_option('-i', '--configinteractive', dest='interactiveConfig',
                      help='generates a default configuration file by interactively ' +
                           'asking the user to input data.',
                      default=False, action='store_true')
    parser.add_option('-q', '--quiet',
                      action='store_false', dest='verbose', default=True,
                      help='don\'t print status messages to stdout')
    parser.add_option('-o', '--overwrite',
                      action='store_true', dest='overwrite', default=False,
                      help='overwrites files. If not present, the pipeline tries ' +
                           'to skip the steps already completed.')

    (options, args) = parser.parse_args()

    if options.interactiveConfig is True:
        options.sampleFile = True

    if len(args) == 0:
        if options.sampleFile is False:
            parser.error('Incorrect number of arguments')
        else:
            PMAP('PMAP.dat', createConfig=True, interactiveConfig=options.interactiveConfig,
                 verbose=options.verbose)
    else:
        configFile = args[0]
        if options.sampleFile is True:
            PMAP(configFile, createConfig=True, interactiveConfig=options.interactiveConfig,
                 verbose=options.verbose)
        else:
            if os.path.exists(configFile):
                PMAP(configFile, createConfig=False, verbose=options.verbose,
                     overwrite=options.overwrite)
            else:
                parser.error('File not found.')
