#!/usr/bin/env python
# encoding: utf-8
"""
getHIFlux.py

Created by Jose Sanchez-Gallego on 29 Aug 2013.
Copyright (c) 2013. All rights reserved.

These functions extract the HI data around FUV regions and
measures the HI flux.
"""

from . import at, pf, pw
from Error import raiseError, raiseWarning
from extractFromCoords import extractFromCoords
import os
import numpy as np
from calcBackground import calcBackground
from . import table
from SExHelpers import call_SEx, read_secat
import glob
import shutil
from astropy import coordinates as coords
from astropy import units


def copySExFiles(logger):
    """
    Copies the default configuration files for SExtractor from the PMAP directory
    to the current path. If the files are already present, it does not overwrite them.
    """

    scriptPath = os.path.dirname(__file__)
    sexFilesPath = scriptPath + '/SExTemplates/*'
    sexFiles = glob.glob(sexFilesPath)
    for sexFile in sexFiles:
        baseName = os.path.basename(sexFile)
        if os.path.exists((baseName)):
            raiseWarning('Sextractor file %s already exists. Not overwriting.' % baseName,
                         logger, newLine=False)
            continue
        shutil.copyfile(sexFile, './' + baseName)

    return


def getHIFlux(configOpts, fluxFUVTable, logger):

    logger.write('Identifying HI patches ...', newLine=True)

    # First, we try to find and read the HI data file
    fluxesHIFile = configOpts['data_hi']
    if os.path.exists(fluxesHIFile) and configOpts['overwrite'] is False:
        try:
            fluxesHI = at.read(fluxesHIFile, delimiter=' ')
            logger.write('... data read from %s' % fluxesHIFile)

            # Very simple and incomplete check for enough data:
            # if len(fluxesHI) != len(fluxFUVTable):
            #     raiseError('Fatal error: number of UV and HI PDRIDs do not match. ' +
            #                'Try regenerating the HI files.')
            # else:
            return fluxesHI
        except:
            raiseError('HI data file %s could not be read.' % fluxesHIFile, logger)

    # If the data file does not exist ...

    imageHI = configOpts['hiImage']
    hduHI = pf.open(imageHI)
    wcsHI = pw.WCS(hduHI[0].header)

    if 'createHIMosaic' not in configOpts:
        createHIMosaic = False
    else:
        createHIMosaic = configOpts['createHIMosaic']

    # We calculate how many pixels are the physical size of the region
    try:
        pixSize = np.array([np.abs(wcsHI.wcs.cd[0, 0]),
                            np.abs(wcsHI.wcs.cd[1, 1])]).mean
    except:
        pixSize = np.array(np.abs(wcsHI.wcs.cdelt[0:2])).mean()

    distance = float(configOpts['distance'])
    physSize = float(configOpts['HIRegSize'])

    scalePc = 2 * distance * np.tan(np.deg2rad(pixSize / 2.)) * 1e6
    size = int(np.ceil(physSize / scalePc))

    dataHI, extractedPaths = extractFromCoords(hduHI, wcsHI, fluxFUVTable, root='HIRegs/HIReg',
                                               extract=True, createMosaic=createHIMosaic,
                                               size=size, overwrite=True, clean=True,
                                               logger=logger)

    # Process HI postage stamp regions.
    # Gets the HI background or automatically calculates it
    if configOpts['hiBackground'] in [None, 'None', 'N', '', '-1', -1.]:
        logger.write('Determining HI background ...', newLine=True)
        hi_bg = calcBackground(imageHI, logger, verbose=configOpts['verbose'])
        configOpts['hiBackground'] = hi_bg
    else:
        hi_bg = configOpts['hiBackground']

    logger.write('Using HI background of %.3e' % hi_bg, newLine=True)

    scale = configOpts['HIscale']  # Do not confuse with scalePc. This one is the
                                   # scale from CPS to cm-2 for the HI map.

    logger.write('Using map scaling of %.3e cm-2 per map unit' % scale)

    # Creates the results table
    fluxesHI = table.Table(None, names=('PDRID', 'HIID', 'RA', 'Dec', 'NHI', 'sNHI'),
                           dtypes=('i4', 'i4', 'S20', 'S20', 'f8', 'f8'))

    logger.write('Copying SExtractor configuration files ...', newLine=True)
    copySExFiles(logger)

    # Warning: up to 999 files
    logger.write('Running SExtractor now ...', newLine=True)

    pdrIDs = sorted(extractedPaths.keys())
    for pdrID in pdrIDs:

        path = extractedPaths[pdrID]

        catFile = call_SEx(path)

        if catFile == 1:
            raiseError('Fatal error. Files %s: SExtractor call failed.' % path, logger)
        else:
            # Moves the SExtractor catalogue just generated to the HI patches directory
            catNewFile = os.path.splitext(path)[0] + '.cat'
            shutil.move(catFile, catNewFile)

            # Reads the catalogue
            patches = read_secat(catNewFile, imageHI)

            # If there has not been a detection
            if patches is None:
                raiseWarning('No suitable patches were found for PDRID %d.' % pdrID, logger, newLine=False)
                # With SExtractor, this is unlikely since it does not pick up the faintest patches
                # fluxesHI.add_row([pdrID, 0, 0, 0, 0])
                continue

            # Patches are not background-subtracted, that will happen now:
            for n in range(len(patches[0])):
                if patches[2][n] > (hi_bg * 2.):
                  # subtract the background and ignore patches that are fainter than hi_bg
                    RA = patches[0][n]
                    Dec = patches[1][n]
                    NHI = (patches[2][n] - hi_bg) * scale
                    # Fix sN_HI to half the background value
                    sNHI = hi_bg * scale * 0.5
                    patchCoords = coords.ICRSCoordinates(ra=RA, dec=Dec, unit=(units.degree, units.degree))
                    RAStr = patchCoords.ra.format(unit=units.hour, precision=2)
                    DecStr = patchCoords.dec.format(unit=units.degree, precision=2)
                    fluxesHI.add_row([pdrID, n, RAStr, DecStr, NHI, sNHI])
                else:
                    raiseWarning('Patch %d for region %d too faint.' % (n, pdrID), logger, newLine=False)

    # Write: HI patches
    at.write(fluxesHI, fluxesHIFile, Writer=at.FixedWidth, delimiter=' ')

    # print "HI patches:"
    # # print "#PDRID, RA, DEC, NHI, sNHI"
    # hifile.write("#PDRID, RA, DEC, NHI, sNHI\n")
    # for n in range(len(hidata)):
    #     # print type(hidata[n].NHI)
    #     # Apparently python 2.6 needs explicit conversion from numpy.float32, so workaround
    #     # fedora bug: changed PDRID:>3 to PDRID:3g
    #     # print "{0.PDRID:3g}, {0.RA:14s}, {0.DEC:14s}, {1:7.5g},
    #     # {2:7.5g}".format(hidata[n], float(hidata[n].NHI), float(hidata[n].sNHI))
    #     hifile.write("{0.PDRID:3g}, {0.RA:14s}, {0.DEC:14s}, {1:7.5g}, {2:7.5g}\n".format(
    #         hidata[n], float(hidata[n].NHI), float(hidata[n].sNHI)))
    #   # hifile.close()
  # End of try/except

    return fluxesHI
