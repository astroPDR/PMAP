#!/usr/bin/env python
# encoding: utf-8
"""
matchFOV.py

Created by José Sánchez-Gallego on 12 Nov 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from .pdrLib import open_image, pf
from astropy.wcs.wcs import FITSFixedWarning
import warnings
import numpy as np
import tempfile

warnings.simplefilter('ignore', FITSFixedWarning)


def getFootPrint(wcs):
    fp = wcs.calcFootprint()
    return np.array([np.min(fp[:, 0]), np.max(fp[:, 0]),
                     np.min(fp[:, 1]), np.max(fp[:, 1])])


def getMatchedArea(wcs1, wcs2):
    """Calculates the common area from two footprints."""

    fp1 = getFootPrint(wcs1)
    fp2 = getFootPrint(wcs2)

    fov = np.zeros(4, dtype=float)

    # Min RA
    if fp1[0] - fp2[0] <= 0.0:
        fov[0] = fp2[0]
    else:
        fov[0] = fp1[0]

    # Max RA
    if fp1[1] - fp2[1] <= 0.0:
        fov[1] = fp1[1]
    else:
        fov[1] = fp2[1]

    # Min dec.
    if fp1[2] - fp2[2] <= 0.0:
        fov[2] = fp2[2]
    else:
        fov[2] = fp1[2]

    # Max dec.
    if fp1[3] - fp2[3] <= 0.0:
        fov[3] = fp1[3]
    else:
        fov[3] = fp2[3]

    fovRect = np.zeros((4, 2), dtype=float)
    fovRect[[0, 1], 0] = fov[1]
    fovRect[[2, 3], 0] = fov[0]
    fovRect[[0, 3], 1] = fov[2]
    fovRect[[1, 2], 1] = fov[3]

    return fovRect


def cutRegion(filename, fov):
    """Cuts an image to a certain FOV and updates the wcs."""

    wcs, data, hdu = open_image(filename)

    coordPix = wcs.wcs_world2pix(fov, 0)

    shape = data.shape

    iMin = np.int(np.min(coordPix[:, 1]))
    iMax = np.int(np.max(coordPix[:, 1]))
    jMin = np.int(np.min(coordPix[:, 0]))
    jMax = np.int(np.max(coordPix[:, 0]))

    if iMin < 0.:
        iMin = 0
    if jMin < 0.:
        jMin = 0
    if iMax > shape[0]:
        iMax = shape[0]
    if jMax > shape[0]:
        jMax = shape[1]

    regionData = data[iMin:iMax, jMin:jMax]
    newHeader = hdu[0].header.copy()

    # We need to update the astrometry of the output image
    if 'CRPIX1' in newHeader.keys():
        newHeader['CRPIX1'] = 1 - ((jMin + 1) - hdu[0].header['CRPIX1'])
        newHeader['CRPIX2'] = 1 - ((iMin + 1) - hdu[0].header['CRPIX2'])
    elif 'LTV1' in newHeader.keys():
        newHeader['LTV1'] = 1 - ((jMin + 1) - hdu[0].header['LTV1'])
        newHeader['LTV2'] = 1 - ((iMin + 1) - hdu[0].header['LTV2'])

    hduPrim = pf.PrimaryHDU(regionData)
    hduPrim.header = newHeader
    hduList = pf.HDUList([hduPrim])

    # tmpFile = tempfile.NamedTemporaryFile(dir='.', suffix='.fits', delete=False)
    # hduPrim.writeto(tmpFile.name)

    return hduList


def maskImage(data, wcs, fov):
    """Mask the region outside the FOV."""

    coordPix = wcs.wcs_world2pix(fov, 0)

    shape = data.shape

    iMin = np.int(np.min(coordPix[:, 1]))
    iMax = np.int(np.max(coordPix[:, 1]))
    jMin = np.int(np.min(coordPix[:, 0]))
    jMax = np.int(np.max(coordPix[:, 0]))

    if iMin < 0.:
        iMin = 0
    if jMin < 0.:
        jMin = 0
    if iMax > shape[0]:
        iMax = shape[0]
    if jMax > shape[0]:
        jMax = shape[1]

    maskedArray = np.ma.array(data, copy=True)

    maskedArray[0:iMin+4, :] = np.ma.masked
    maskedArray[iMax-4:, :] = np.ma.masked
    maskedArray[:, 0:jMin+4] = np.ma.masked
    maskedArray[:, jMax-4:] = np.ma.masked

    return maskedArray


def matchFOV(image1, image2):
    """Calculates the matching footprint of a pair of images.

    Calculates the footprint of both images and the ommon area. Returns an
    array with the corners of the common area (similar to astropy.wcs.WCS.calcFootprint).

    """

    wcs1, data1, hdu1 = open_image(image1)
    wcs2, data2, hdu2 = open_image(image2)

    fovRect = getMatchedArea(wcs1, wcs2)

    return fovRect
