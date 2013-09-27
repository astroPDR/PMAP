#!/usr/bin/env python
# encoding: utf-8
"""
extractRegions.py

Created by José Ramón Sánchez-Gallego on 2011-03-28.
Copyright (c) 2011. All rights reserved.

This script extracts the images
Needs: pyfits, numpy
Usage: extractRegions.py [-m|--mosaic] [-e|--extend] image mask output
        image: Original FITS image of the object
        mask: The file with the masked regions
        output: The root of the output extracted regions
        -e|--extend: Number of pixels the extracted region will be extended
        -m|--mosaic: creates a mosaic with the extracted regions

"""

import numpy as np
import os
from . import pf
from Regions import *

def flush(text):
    import sys
    sys.stdout.write(text)
    sys.stdout.flush()
    return


def run(cmd):
    import subprocess
    subprocess.call(cmd, shell=True, stdout=open('/dev/null', 'w'), stderr=subprocess.STDOUT)
    return


def extractRegions(image, mask, output, extend=2, mosaic=False, mosaicBack=0):

# Loads the original image and mask
    hduImage = pf.open(image)
    hduMask = pf.open(mask)

# Loads the regions
    flush('Loading regions ... ')
    regions = RegionSet(hduMask[0].data, image=hduImage[0].data)
    print 'done'

# For each region, calculates its limits, trims the original image and
# saves the new file

    def saveData(data, outFormat, message):

        nRegions = regions.nRegs
        nZeros = int(np.floor(np.log10(nRegions) + 1))

        for ii,idReg in enumerate(regions.idRegions):

            iMin = np.int(np.min(regions.Regions[idReg].pixels[:,0])-extend)
            # iMax = np.int(np.max(regions.Regions[idReg].pixels[:,0])+extend)
            jMin = np.int(np.min(regions.Regions[idReg].pixels[:,1])-extend)
            # jMax = np.int(np.max(regions.Regions[idReg].pixels[:,1])+extend)

            region = regions.Regions[idReg].getSlice(data, extend=extend)[0]
            hduPrim = pf.PrimaryHDU(region)

# We need to update the astrometry of the output image
            newHeader = hduImage[0].header.copy()
            if 'CRPIX1' in newHeader.keys():
                newHeader['CRPIX1'] = 1 - ((jMin + 1) - hduImage[0].header['CRPIX1'])
                newHeader['CRPIX2'] = 1 - ((iMin + 1) - hduImage[0].header['CRPIX2'])
            elif 'LTV1' in newHeader.keys():
                newHeader['LTV1'] = 1 - ((jMin + 1) - hduImage[0].header['LTV1'])
                newHeader['LTV2'] = 1 - ((iMin + 1) - hduImage[0].header['LTV2'])

            hduPrim.header = newHeader
            hduReg = pf.HDUList([hduPrim])

            outName = outFormat.format(output, str(idReg).zfill(nZeros))
            if os.path.exists(outName): os.remove(outName)
            if not os.path.exists('{0}_Regions'.format(output)): os.mkdir('{0}_Regions'.format(output))
            hduReg.writeto(outName, output_verify='ignore')

            del hduPrim
            del hduReg

            percent = np.int(np.rint((ii) * 100. / nRegions))
            flush('\r%s ... %d%% ' % (message, percent))

    flush('\rSaving image regions ... 0%')
    saveData(regions.Image, '{0}_Regions/{0}_Reg{1}.fits', 'Saving image regions')
    flush('\rSaving image regions ... Done\n')
    flush('\rSaving mask regions ... 0%')
    saveData(regions.DataMask, '{0}_Regions/{0}_Msk{1}.fits', 'Saving mask regions')
    flush('\rSaving mask regions ... Done\n')


    if mosaic:

        regSize = np.array([regions.Regions[id].sliceMask.shape for id in regions.idRegions]) - 2
        maxSizes = np.max(regSize, axis=0)
        if regions.nRegs <= 8:
            nCols = 8
            nRows = 1
        else:
            nCols = int(np.sqrt(regions.nRegs)) + 1
            if nCols**2 >= regions.nRegs:
                nRows = nCols
            else:
                nRows = nCols + 1

        def saveMosaic(data, outName):
            frame = np.zeros((maxSizes+4*extend) * np.asarray([nRows, nCols])) + mosaicBack

            ii = 1
            jj = 0
            for nn,idReg in enumerate(regions.idRegions):

                iiStart = frame.shape[0] - ((maxSizes[0] + 4*extend) * ii) + extend
                jjStart = ((maxSizes[1] + 4*extend) * jj) + extend

                slice = regions.Regions[idReg].getSlice(data, extend=extend)[0]
                frame[iiStart:iiStart+slice.shape[0], jjStart:jjStart+slice.shape[1]] = slice

                jj += 1
                if jj > nCols-1:
                    jj = 0
                    ii += 1

            hdu = pf.PrimaryHDU(frame)
            hduList = pf.HDUList([hdu])
            if os.path.exists(outName): os.remove(outName)
            hduList.writeto(outName, output_verify='ignore')

        flush('Creating and saving mosaics ... ')
        saveMosaic(regions.Image, output + '_Mosaic.fits')
        saveMosaic(regions.DataMask, output + '_Mosaic_Msk.fits')
        print 'done'

    return


# Parses the command line options and call the main extractRegions routine
if __name__ == '__main__':

    from optparse import OptionParser

    parser = OptionParser()

    print

    usage = 'usage: %prog [options] image mask output'
    parser = OptionParser(usage=usage)
    parser.add_option('-x', '--extend', dest='extend',
        help='number of pixels the extracted region will be extended', type=int, default=2)
    parser.add_option('-m', '--mosaic', dest='mosaic', action='store_true', default=False,
        help='creates a mosaic with the extracted regions')

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error('incorrect number of arguments\n')
    elif len(args) == 1:
        image = args[0]
        mask = os.path.splitext(image)[0] + '_Msk.fits'
        print 'Assuming Mask file: {0}'.format(mask)
        output = os.path.splitext(image)[0]
        print 'Assuming Output file: {0}\n'.format(output)
    elif len(args) == 2:
        image = args[0]
        mask = args[1]
        output = os.path.splitext(image)[0]
        print 'Assuming Output: {0}\n'.format(output)
    elif len(args) == 3:
        image = args[0]
        mask = args[1]
        output = args[2]

    if False in map(os.path.exists, [image, mask]):
        parser.error('Either the image or the mask does not exist\n')

    extractRegions(image, mask, output, extend=options.extend, mosaic=options.mosaic)
    print

