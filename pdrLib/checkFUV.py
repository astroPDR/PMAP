#!/usr/bin/env python
# encoding: utf-8
"""
checkFUV.py

Created by José Ramón Sánchez-Gallego on 2011-03-08.
Copyright (c) 2011. All rights reserved.

This script localized regions that are "touching" and checks them to make sure
that they are two different regions. To do that:
(1) Finds the peak position for each test region
(2) Localizes the pixels joining the peak positions with a straight line
(3) Calculates the intensity profile across that line
(4) If certain conditions are not fulfilled, both test regions are joined
in a single region

Modified on 2013-04-16 by José Sánchez-Gallego
  -- General cleanup. No major changes.

"""

import numpy as np
import os
from pdrLib import pf
import Regions
from Error import raiseError


def calcAngDist(coords1, coords2):
# Taken from astCoords

    RADeg1 = coords1[0]
    RADeg2 = coords2[0]
    decDeg1 = coords1[1]
    decDeg2 = coords2[1]

    cRA = np.radians(RADeg1)
    cDec = np.radians(decDeg1)

    gRA = np.radians(RADeg2)
    gDec = np.radians(decDeg2)

    cosC = (np.sin(gDec) * np.sin(cDec)) + (np.cos(gDec) * np.cos(cDec) * np.cos(gRA - cRA))
    xx = (np.cos(cDec)*np.sin(gRA-cRA))/cosC
    yy = ((np.cos(gDec) * np.sin(cDec)) - (np.sin(gDec) * np.cos(cDec) * np.cos(gRA - cRA))) / cosC
    rr = np.degrees(np.sqrt(xx * xx + yy * yy))

    return rr


def checkFUV(options, logger):

    fuvImage = options['fuvImage']
    fuvMask = options['fuvMaskFile']

    logger.write('Checking FUV regions ... ', newLine=True, doLog=False)

    logger.write('Gathering regions ... ', doLog=False)

    hduMask = pf.open(fuvMask)
    hduImage = pf.open(fuvImage)
    # nRegs = np.max(hduMask[0].data)
    regions = Regions.RegionSet(hduMask[0].data, image=hduImage[0].data)

    logger.write('Getting adjacent regions ... ', doLog=False)
    adjacentRegions = regions.getAdjacentRegs()

    logger.write('Processing adjacent regions ... ', doLog=False)

    while adjacentRegions != []:

        pair = adjacentRegions[0]
        fluxA = regions.Regions[pair[0]].getFlux()
        fluxB = regions.Regions[pair[1]].getFlux()
        peakA = regions.Regions[pair[0]].peak
        peakB = regions.Regions[pair[1]].peak

        if options['fluxContrast'] is not False:
            contrast = float(np.max([fluxA, fluxB]) / np.min([fluxA, fluxB]))
            minRegion = pair[np.argmin([fluxA, fluxB])]
            maxRegion = pair[np.argmax([fluxA, fluxB])]
            if options['fluxContrast'] is True or contrast > options['fluxContrast']:
                if options['joinRegions']:
                    regions.joinRegions(maxRegion, minRegion)
                    logger.write('Adjacent regions %d and %d joined. Flux contrast=%.1f' %
                                 (maxRegion, minRegion, contrast), doPrint=False)
                else:
                    regions.deleteReg(minRegion)
                    logger.write('Adjacent regions %d deleted. Flux contrast=%.1f' %
                                 (minRegion, contrast), doPrint=False)
                adjacentRegions = regions.getAdjacentRegs()
                continue

        if options['peakContrast'] is not False:
            contrast = float(np.max([peakA, peakB]) / np.min([peakA, peakB]))
            minRegion = pair[np.argmin([peakA, peakB])]
            maxRegion = pair[np.argmax([peakA, peakB])]
            if options['peakContrast'] is True or contrast > options['peakContrast']:
                if options['joinRegions']:
                    regions.joinRegions(maxRegion, minRegion)
                    logger.write('Adjacent regions %d and %d joined. Peak contrast=%.1f' %
                                 (maxRegion,  minRegion, contrast), doPrint=False)
                else:
                    regions.deleteReg(minRegion)
                    logger.write('Adjacent regions %d deleted. Peak contrast=%.1f' %
                                 (minRegion, contrast), doPrint=False)
                adjacentRegions = regions.getAdjacentRegs()
                continue

        del adjacentRegions[0]

    logger.write('Getting close regions ... ', newLine=True, doLog=False)

    if options['scale'] is not None:
        scalePc = options['scale']
    else:
        try:
            try:
                scaleDeg = np.abs(hduImage[0].header['CD1_1'])
            except:
                scaleDeg = np.abs(hduImage[0].header['CDELT1'])
        except:
            raiseError('Pixel scale cannot be calculated. Use the scale parameter.')

        try:
            scalePc = 2.0 * options['distance'] * 1e6 * np.tan(0.5 * scaleDeg * np.pi / 180.)
        except:
            raiseError('Distance to the galaxy not defined.')

    minNumPixels = options['minDistance'] / scalePc

    closeRegions = regions.getCloseRegs(minNumPixels, adjacent=False)

    logger.write('Processing close regions ... ', newLine=True, doLog=False)

    while closeRegions != []:

        pair = closeRegions[0]

        fluxA = regions.Regions[pair[0]].getFlux()
        fluxB = regions.Regions[pair[1]].getFlux()
        peakA = regions.Regions[pair[0]].peak
        peakB = regions.Regions[pair[1]].peak

        if options['alwaysRemoveClose']:
            minRegion = pair[np.argmin([fluxA, fluxB])]
            regions.deleteReg(minRegion)
            closeRegions = regions.getCloseRegs(minNumPixels, adjacent=False)
            logger.write('Close region %d removed' % minRegion, doPrint=False)
            continue

        if options['fluxContrast'] is not False:
            contrast = float(np.max([fluxA, fluxB]) / np.min([fluxA, fluxB]))
            minRegion = pair[np.argmin([fluxA, fluxB])]
            maxRegion = pair[np.argmax([fluxA, fluxB])]
            if (options['fluxContrast'] is True) or (contrast > options['fluxContrast']):
                regions.deleteReg(minRegion)
                closeRegions = regions.getCloseRegs(minNumPixels, adjacent=False)
                logger.write('Close region %s removed. Flux contrast=%.1f' %
                             (minRegion, contrast), doPrint=False)
                continue

        if options['peakContrast'] is not False:
            contrast = float(np.max([peakA, peakB]) / np.min([peakA, peakB]))
            minRegion = pair[np.argmin([peakA, peakB])]
            maxRegion = pair[np.argmax([peakA, peakB])]
            if (options['peakContrast'] is True) or (contrast > options['peakContrast']):
                regions.deleteReg(minRegion)
                closeRegions = regions.getCloseRegs(minNumPixels, adjacent=False)
                logger.write('Close region %d removed. Peak contrast=%.1f' %
                             (minRegion, contrast), doPrint=False)
                continue

        del closeRegions[0]

    logger.write('Saving new mask ... ', newLine=True, doLog=False)

    hdu = pf.PrimaryHDU(regions.DataMask)
    hdu.header = hduMask[0].header
    hduList = pf.HDUList([hdu])
    if os.path.exists(options['fuvMaskFileRej']):
        os.remove(options['fuvMaskFileRej'])
    hduList.writeto(options['fuvMaskFileRej'], output_verify='ignore')

    return
