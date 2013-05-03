#!/usr/bin/env python
# encoding: utf-8
"""
extractFromCoords.py

Created by José Ramón Sánchez-Gallego on 2011-05-19.
Copyright (c) 2011. All rights reserved.

This scripts accepts a coordinates file in DS9 format and extracts
regions from a FITS file.

"""

import sys, os
import numpy as np
# import pyregion
import pywcs, pyfits

def extractFromCoords(image, ds9File, output, imageDS9=None, size=20, mosaic=False, mosaicBack=0):

    sys.stdout.write('Reading DS9 file ... ')
    sys.stdout.flush()

    try:
        regions = pyregion.parse(open(ds9File, 'r').read())
    except:
        print
        return

    print 'done'

    sys.stdout.write('Loading image ... ')
    sys.stdout.flush()

    try:
        hdu = pyfits.open(image)
    except:
        print 'failed'
        return

    try:
        wcs = pywcs.WCS(hdu[0].header)
    except:
        print 'failed. No WCS for the main image'
        print
        return

    if imageDS9 != None:
        try:
            wcsDS9 = pywcs.WCS(pyfits.getheader(imageDS9))
        except:
            wcsDS9 = None

    print 'done'


    sys.stdout.write('Reading coordinates ... ')
    sys.stdout.flush()

    coordsSky = []
    for region in regions:
        system = region.coord_format
        coords = region.coord_list[0:2]
        if system == 'image':
            if wcsDS9 == None:
                print 'failed. No WCS for the reference image.'
                print
                return
            else:
                coords = wcsDS9.wcs_pix2sky(np.array([coords])-1.,0)[0]

        coordsSky.append(coords)

    coordsSky = np.array(coordsSky)
    print 'done'

    sys.stdout.write('Saving regions ... ')
    sys.stdout.flush()

    nRegions = coordsSky.shape[0]
    nZeros = 3 #int(np.floor(np.log10(nRegions) + 1)) #flexible is nice but code later on is not

    mosaicData = []
    for ii,coordSky in enumerate(coordsSky):

        coordPix = wcs.wcs_sky2pix(np.array([coordSky]),0)[0] + 1

        iMin = np.int(coordPix[1] - size)
        iMax = np.int(coordPix[1] + size)
        jMin = np.int(coordPix[0] - size)
        jMax = np.int(coordPix[0] + size)

        regionData = hdu[0].data[iMin:iMax, jMin:jMax]
        hduPrim = pyfits.PrimaryHDU(regionData)
        newHeader = hdu[0].header.copy()

# We need to update the astrometry of the output image
        if wcs != None:
            if 'CRPIX1' in newHeader.keys():
                newHeader['CRPIX1'] = 1 - ((jMin + 1) - hdu[0].header['CRPIX1'])
                newHeader['CRPIX2'] = 1 - ((iMin + 1) - hdu[0].header['CRPIX2'])
            elif 'LTV1' in newHeader.keys():
                newHeader['LTV1'] = 1 - ((jMin + 1) - hdu[0].header['LTV1'])
                newHeader['LTV2'] = 1 - ((iMin + 1) - hdu[0].header['LTV2'])

        hduPrim.header = newHeader
        hduReg = pyfits.HDUList([hduPrim])

        outName = '{0}_Regions/{0}_Reg{1}.fits'.format(output, str(ii+1).zfill(nZeros))
        if os.path.exists(outName): os.remove(outName)
        if not os.path.exists('{0}_Regions'.format(output)):
            os.mkdir('{0}_Regions'.format(output))
        hduReg.writeto(outName, output_verify='ignore')

        del hduPrim
        del hduReg

        percent = np.int(np.rint((ii) * 100. / nRegions))
        sys.stdout.write('\rSaving regions ... %d%% ' % percent)
        sys.stdout.flush()

        if mosaic == True: mosaicData.append(regionData)

    print '\rSaving regions ... done'


    if mosaic:

        maxSizes = 2.0*size
        if nRegions <= 8:
            nCols = 8
            nRows = 1
        else:
            nCols = int(np.sqrt(nRegions)) + 1
            if nCols**2 >= nRegions:
                nRows = nCols
            else:
                nRows = nCols + 1

        def saveMosaic(data, outName):
            extend = 2.
            frame = np.zeros((maxSizes+4*extend) * np.asarray([nRows, nCols])) + mosaicBack

            ii = 1
            jj = 0
            for nn,regData in enumerate(data):

                iiStart = frame.shape[0] - ((maxSizes + 4*extend) * ii) + extend
                jjStart = ((maxSizes + 4*extend) * jj) + extend

                frame[iiStart:iiStart+maxSizes, jjStart:jjStart+maxSizes] = regData

                jj += 1
                if jj > nCols-1:
                    jj = 0
                    ii += 1

            hduMosaic = pyfits.PrimaryHDU(frame)
            hduListMosaic = pyfits.HDUList([hduMosaic])
            if os.path.exists(outName): os.remove(outName)
            hduListMosaic.writeto(outName, output_verify='ignore')

        sys.stdout.write('Creating and saving mosaics ... ')
        sys.stdout.flush()
        saveMosaic(mosaicData, output + '_Mosaic.fits')
        print 'done'


    print

    return


if __name__ == '__main__':

    from optparse import OptionParser

    parser = OptionParser()

    print

    usage = 'usage: %prog [options] image ds9File output'
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--size', dest='size',
        help='width, in pixels, of the squared regions to extract', type=int, default=20)
    parser.add_option('-d', '--ds9', dest='imageDS9',
            help='image to which the ds9File is referred, if the coordinate system is not fk5',
            type=str, default=None)
    parser.add_option('-m', '--mosaic', dest='mosaic', action='store_true', default=False,
        help='creates a mosaic with the extracted regions')

    (options, args) = parser.parse_args()

    if len(args) <= 1:
        parser.error('incorrect number of arguments\n')
    elif len(args) == 2:
        image = args[0]
        ds9File = args[1]
        output = os.path.splitext(image)[0]
        print 'Assuming Output: {0}\n'.format(output)
    elif len(args) == 3:
        image = args[0]
        ds9File = args[1]
        output = args[2]

    if False in map(os.path.exists, [image, ds9File]):
        parser.error('Either the image or the mask does not exist\n')

    extractFromCoords(image, ds9File, output, size=options.size,
        mosaic=options.mosaic, imageDS9=options.imageDS9)
