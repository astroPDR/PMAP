#!/usr/bin/env python
# encoding: utf-8
"""
getRegions.py

Created by José Ramón Sánchez-Gallego on 2011-01-17.
Modified by José Ramón Sánchez-Gallego on 2011-07-04.
JSH 2011-10-12: changed np.max and np.min to np.nanmax and np.nanmin respectively

Modified by José Ramón Sánchez-Gallego on 2013-04-13.
    -- Clean-up and small bug fixes
    -- Now it uses the native Python CLFind.

Copyright (c) 2011. All rights reserved.

This file contains several routines to read the configuration
file and to perform the selection of regions for the FUV and HI
images.

"""

import os
import numpy as np
from Error import raiseError
from calcBackground import calcBackground
from clfind2d import clfind2d
from checkFUV import checkFUV
from createCatalog import createCatalog
from pdrLib import open_image
from matchFOV import maskImage


def checkCLFindOptions(root, configOpts):
    everyThingOK = True
    try:
        if configOpts[root + 'Background'] is None:
            everyThingOK = False
        if configOpts[root + 'hLevel'] is None:
            everyThingOK = False
        if configOpts[root + 'nLevels'] is None:
            everyThingOK = False
        if configOpts[root + 'Sigma'] is None:
            everyThingOK = False
        if configOpts[root + 'MinPixels'] is None:
            everyThingOK = False
    except:
        everyThingOK = False
    return everyThingOK


def calcParameters(root, configOpts, logger, sigmaSep=3., fov=None):

    image = configOpts[root + 'Image']
    sigma = configOpts[root + 'Sigma']

    if configOpts[root + 'Background'] in [None, 'None', 'N', '', '-1', -1.]:
        logger.write('Determining FUV background ...', newLine=True)
        background = calcBackground(image, logger,
                                    verbose=configOpts['verbose'], fov=fov)
        configOpts[root + 'Background'] = background
    else:
        background = configOpts[root + 'Background']

    lLevel = background * sigma
    configOpts[root + 'lLevel'] = lLevel

    if configOpts[root + 'hLevel'] in [None, 'None', 'N', '', '-1', -1.]:

        wcs, data, hdu = open_image(image)
        if fov is not None:
            data = maskImage(data, wcs, fov)

        #hLevel = np.nanmax(data) - background
        hLevel = np.max(data) - background

        configOpts[root + 'hLevel'] = hLevel

    hLevel = configOpts[root + 'hLevel']

    if configOpts[root + 'nLevels'] in [None, 'None', 'N', '', '-1', -1.]:
        nLevels = np.floor((hLevel - lLevel)/(sigmaSep * background))
        configOpts[root + 'nLevels'] = nLevels
        configOpts[root + 'Levels'] = np.arange(lLevel, hLevel, sigmaSep * background)
    else:
        nLevels = configOpts[root + 'nLevels']
        configOpts[root + 'Levels'] = np.arange(lLevel, hLevel, (hLevel-lLevel)/nLevels)

    return


def printCLFindOptions(root, configOpts, logger):

    logger.write('Running CLFind with the following parameters: ', newLine=True)

    for param in ['Image', 'Background', 'Sigma', 'lLevel', 'hLevel', 'nLevels']:
        logger.write(root + param + ' = ' + str(configOpts[root + param]))

    return


def getYN(text, returnOpt=None):
    opt = 'x'
    validKeys = {'y': True, 'n': False}
    if returnOpt is not None:
        validKeys[''] = returnOpt
    while opt.lower() not in validKeys.keys():
        opt = raw_input(text)
    return validKeys[opt]


def getFloat(text, default):
    value = None
    while type(value) not in [float, int]:
        value = raw_input(text)
        try:
            if value == '':
                return default
            value = float(value)
        except:
            pass
    return value


def getRegions(configOpts, logger, fov=None):

    #####################################################################
    ##########################  FUV processing ##########################

    fuvImage = configOpts['fuvImage']

    logger.write('FUV Image: %s' % fuvImage, newLine=True)

    if not os.path.exists(fuvImage):
        raiseError('FUV image does not exist.', logger)

    if (configOpts['fuvMaskFile'].lower() == 'default') or (configOpts['fuvMaskFile'] is None):
        configOpts['fuvMaskFile'] = os.path.splitext(configOpts['fuvImage'])[0] + '_Mask.fits'

    if fov is not None:
        maskFile = os.path.realpath(configOpts['fuvMaskFile'])
        trimmedFUVImage = os.path.dirname(maskFile) + '/Trimmed_' + \
            os.path.basename(maskFile)
        if os.path.exists(trimmedFUVImage):
            configOpts['fuvImage'] = trimmedFUVImage
            logger.write('FUV trimmed image is {0}'.format(os.path.basename(trimmedFUVImage)),
                         newLine=True, doPrint=True)

    doCLFindFUV = True
    if os.path.exists(configOpts['fuvMaskFile']):
        if configOpts['overwrite'] is True:
            doCLFindFUV = True
        else:
            doCLFindFUV = getYN('\nFound %s. Rerun CLFind? [y/N] ' %
                                os.path.basename(configOpts['fuvMaskFile']),
                                returnOpt=False)

    if doCLFindFUV is False:
        logger.write('FUV mask not created. Using previously generated mask %s.' %
                     os.path.basename(configOpts['fuvMaskFile']))

    else:

        configOpts['overwrite'] = True  # From this point on we have to redo everything

        if configOpts['fuvInteractive'] is False:
            logger.write('FUV parameters loaded')

        else:

            logger.write('FUV: Running in interactive mode ... ', newLine=True)

            wcs, data, hdu = open_image(fuvImage)
            if fov is not None:
                data = maskImage(data, wcs, fov)

            logger.write('Image statistics', newLine=True)
            logger.write('Max: %9.3e' % np.nanmax(data))
            logger.write('Min: %9.3e' % np.nanmin(data))
            logger.write('Mean: %9.3e' % np.mean(data))
            logger.write('Median: %9.3e' % np.median(data))

            configOpts['fuvBackground'] = getFloat('Background level [-1 for automatic]: ', -1)
            configOpts['fuvSigma'] = getFloat('Sigma [3]: ', 3)
            configOpts['fuvhLevel'] = getFloat('Highest level [-1 for automatic]: ', -1)
            configOpts['fuvnLevels'] = getFloat('Number of levels  [-1 for automatic]: ', -1)
            configOpts['fuvMinPixels'] = getFloat('Minimum number of pixels [20]: ', 20)

        calcParameters('fuv', configOpts, logger, fov=fov)
        printCLFindOptions('fuv', configOpts, logger)

        resultCheck = checkCLFindOptions('fuv', configOpts)

        if resultCheck is False:
            raiseError('Some of the parameters are incorrect. ' +
                       'Please, review the configuration file.',
                       logger)

        clfind2d(fuvImage, configOpts['fuvMaskFile'], configOpts['fuvLevels'],
                 nPixMin=configOpts['fuvMinPixels'], verbose=configOpts['verbose'],
                 extend=fov)

        logger.write('CLFind log can be found in %s.log' %
                     configOpts['fuvMaskFile'], newLine=True, doPrint=False)

    #####################################################################
    ##########################  FUV rejection ###########################

    # We calculate how many pixels are the physical size of the region
    physSize = float(configOpts['HIRegSize'])
    scalePc = configOpts['HI_PcPerPixel']
    minStampSize = int(np.ceil(physSize / scalePc))

    if configOpts['checkFUV']:
        if (configOpts['fuvMaskFileRej'] is None) or \
                (configOpts['fuvMaskFileRej'].lower() == 'default'):
            configOpts['fuvMaskFileRej'] = os.path.splitext(configOpts['fuvImage'])[0] + \
                '_RejMsk.fits'

        doRejection = True
        if os.path.exists(configOpts['fuvMaskFileRej']) and doCLFindFUV is False:
            if configOpts['overwrite'] is False:
                doRejection = getYN('\nFound %s. Redo FUV rejection? [y/N] ' %
                                    os.path.basename(configOpts['fuvMaskFileRej']),
                                    returnOpt=False)

        if doRejection is True:
            configOpts['overwrite'] = True  # From this point on we have to redo everything
            logger.write('FUV rejection', newLine=True)
            checkFUV(configOpts, logger, minStampSize=minStampSize)
            logger.write('New mask saved as %s' % configOpts['fuvMaskFileRej'],
                         newLine=True)
        else:
            logger.write('Not doing FUV rejection. Using image %s' %
                         configOpts['fuvMaskFileRej'])

    #####################################################################
    ######################  FUV regions catalogue #######################

    logger.write('Creating FUV catalogs ... ', newLine=True)

    if configOpts['checkFUV'] is True:
        fuvMask = configOpts['fuvMaskFileRej']
    else:
        fuvMask = configOpts['fuvMaskFile']

    configOpts['finalFUVMask'] = fuvMask  # Saves the final FUV mask to the config

    votRegsFile = os.path.splitext(fuvMask)[0] + '.vot'
    ds9RegsFile = os.path.splitext(fuvMask)[0] + '.reg'
    peaksFile = os.path.splitext(fuvMask)[0] + '_Peaks.dat'
    configOpts['fuvPeaksFile'] = peaksFile

    if False not in map(os.path.exists, [votRegsFile, ds9RegsFile, peaksFile]) and \
            configOpts['overwrite'] is False:
        logger.write('FUV catalogues already exist.', newLine=True)
        if os.path.exists(votRegsFile):
            logger.write('VOTable catalogue: %s' % votRegsFile)
        logger.write('DS9 catalogue: %s' % ds9RegsFile)
        logger.write('Peaks file: %s' % peaksFile)
    else:
        createCatalog(fuvImage, fuvMask, votRegsFile, ds9RegsFile, peaksFile,
                      logger, ellipse=False)

    return
