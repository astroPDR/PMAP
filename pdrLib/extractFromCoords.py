#!/usr/bin/env python
# encoding: utf-8
"""
extractFromCoords.py

Created by José Ramón Sánchez-Gallego on 2011-05-19.
Copyright (c) 2011. All rights reserved.

Modified 2013-08-28 by José R. Sánchez-Gallego.

This scripts accepts a pyfits HDU and FUV flux table and
extracts the regions around the FUV peaks.

"""

import os
import numpy as np
from astropy import coordinates as coord
from . import pf
from Error import raiseWarning


def extractFromCoords(hduHI, wcsHI, fluxTable, root='HIRegs/HIReg',
                      extract=False, createMosaic=False, mosaicBack=0,
                      size=20, overwrite=False, logger=None, clean=False,
                      astropyVersion='0.2.4'):

    # If clean is True, we delete the path before doing anything
    dirName = os.path.dirname(root)
    if clean is True:
        if os.path.exists(dirName):
            import shutil
            shutil.rmtree(dirName)

    outPutData = {}
    extractedPaths = {}
    mosaicData = []
    for nn, peak in enumerate(fluxTable):

        ID = int(peak['PDRID'])
        RA = peak['RA']
        Dec = peak['Dec']
        peakCoords = coord.ICRSCoordinates('%s %s' % (RA, Dec))

        if astropyVersion == '0.2.4':
            coordPix = wcsHI.wcs_world2pix(np.array([[peakCoords.ra.degrees,
                                                      peakCoords.dec.degrees]]), 0)[0] + 1
        else:
            coordPix = wcsHI.wcs_world2pix(np.array([[peakCoords.ra.degree,
                                                      peakCoords.dec.degree]]), 0)[0] + 1

        iMin = np.int(coordPix[1] - size)
        iMax = np.int(coordPix[1] + size)
        jMin = np.int(coordPix[0] - size)
        jMax = np.int(coordPix[0] + size)

        regionData = hduHI[0].data[iMin:iMax, jMin:jMax]
        outPutData[ID] = regionData
        if createMosaic is True:
            mosaicData.append(regionData)

        # If extract is True, we have to save regionData as a new
        # FITS image

        # First, we check if the image exists and if it must be overwritten
        if extract is True:
            path = root + '%04d.fits' % ID
            if os.path.exists(path):
                if overwrite is False:
                    if nn == 0:
                        raiseWarning('Image %s already exists. Not overwriting.' % path,
                                     logger, doPrint=True, doLog=True, newLine=False)
                        raiseWarning('The remaining warning of this type will be hidden.',
                                     logger, doPrint=True, doLog=False, newLine=False)
                    else:
                        raiseWarning('Image %s already exists. Not overwriting.' % path,
                                     logger, doPrint=False, doLog=True, newLine=False)
                    extractedPaths[ID] = path
                    continue
                else:
                    os.remove(path)

            if not os.path.exists(dirName):
                os.mkdir(dirName)

            # Now, we create the new FITS image
            hduPrim = pf.PrimaryHDU(regionData)
            newHeader = hduHI[0].header.copy()

            # We need to update the astrometry of the output image
            if wcsHI is not None:
                if 'CRPIX1' in newHeader.keys():
                    newHeader['CRPIX1'] = 1 - ((jMin + 1) - hduHI[0].header['CRPIX1'])
                    newHeader['CRPIX2'] = 1 - ((iMin + 1) - hduHI[0].header['CRPIX2'])
                elif 'LTV1' in newHeader.keys():
                    newHeader['LTV1'] = 1 - ((jMin + 1) - hduHI[0].header['LTV1'])
                    newHeader['LTV2'] = 1 - ((iMin + 1) - hduHI[0].header['LTV2'])

            hduPrim.header = newHeader
            hduReg = pf.HDUList([hduPrim])

            hduReg.writeto(path, output_verify='ignore')
            extractedPaths[ID] = path

            del hduPrim
            del hduReg

    if createMosaic is True:

        mosaicPath = os.path.dirname(root) + '/Mosaic.fits'
        if os.path.exists(mosaicPath):
            if overwrite is False:
                raiseWarning('Mosaic %s already exists. Not overwriting.' % mosaicPath,
                             logger, doPrint=True, doLog=True, newLine=False)
                return outPutData, extractedPaths

        maxSizes = 2.0 * size
        nRegions = len(fluxTable)
        if nRegions <= 8:
            nCols = 8
            nRows = 1
        else:
            nCols = int(np.sqrt(nRegions)) + 1
            if nCols ** 2 >= nRegions:
                nRows = nCols
            else:
                nRows = nCols + 1

        def saveMosaic(data, outName):
            extend = 2.
            frame = np.zeros((maxSizes + 4 * extend) * np.asarray([nRows, nCols])) + mosaicBack

            ii = 1
            jj = 0
            for nn, regData in enumerate(data):

                iiStart = frame.shape[0] - ((maxSizes + 4 * extend) * ii) + extend
                jjStart = ((maxSizes + 4 * extend) * jj) + extend

                frame[iiStart:iiStart + maxSizes, jjStart:jjStart + maxSizes] = regData

                jj += 1
                if jj > nCols - 1:
                    jj = 0
                    ii += 1

            hduMosaic = pf.PrimaryHDU(frame)
            hduListMosaic = pf.HDUList([hduMosaic])
            if os.path.exists(outName):
                    os.remove(outName)
            hduListMosaic.writeto(outName, output_verify='ignore')

        saveMosaic(mosaicData, mosaicPath)

    return outPutData, extractedPaths
