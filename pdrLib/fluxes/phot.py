#!/usr/bin/env python
# encoding: utf-8
"""
phot.py

Created by José Sánchez-Gallego on 3 May 2013.
Copyright (c) 2013. All rights reserved.

"""

import numpy as np
import astropy.io.ascii as at


def rebin(aa, shape):
    """
    Resizes an array to a certain shape
    """
    sh = shape[0], aa.shape[0] // shape[0], shape[1], aa.shape[1] // shape[1]
    return aa.reshape(sh).mean(-1).mean(1)


def calcMask(data, xCentre, yCentre, r1, r2=None, refine=False, nRefine=10):
    """
    Returns a mask with the pixels included in the aperture. If r2 is defined,
    the apertures is annular. If refine, the data is re-meshed by a factor
    nRefine.
    """

    if r2 is None:
        rMax = int(r1)
    else:
        rMax = int(np.max([r1, r2]))

    # Grid of coordinates
    if refine is True:
        step = 1. / nRefine
    else:
        step = 1.

    # Calculates the coordinate grid
    ii = yCentre + np.arange(-rMax-2, rMax+2, step)
    jj = xCentre + np.arange(-rMax-2, rMax+2, step)
    coords = np.dstack(np.meshgrid(ii, jj)).reshape(-1, 2)

    # Calculates distances
    distances = (coords[:, 0]-yCentre)**2 + (coords[:, 1]-xCentre)**2
    distances = distances.reshape(ii.size, jj.size)

    mask = np.zeros(distances.shape)

    # Determines pixels within aperture and set them to 1.
    if r2 is None:
        mask[distances <= (r1**2)] = 1
    else:
        mask[(distances >= (r1**2)) & (distances <= (r2**2))] = 1

    # From the data, selects the area matching the mask
    dataChunk = data[np.round(yCentre-1)-rMax-2:np.round(yCentre-1)+rMax+2,
                     np.round(xCentre-1)-rMax-2:np.round(xCentre-1)+rMax+2]

    # Rebins the mask to its original size
    mask = rebin(mask, dataChunk.shape)
    # Pixels outside aperture are set to nan
    mask[mask == 0] = np.nan

    del distances, coords, ii, jj

    return dataChunk, mask


def getRejPixels(origData, nSigma=3., nMaxIter=50):
    """
    Rejects pixels outside a certain range of sigma
    """

    data = np.array(origData, copy=True)
    nSrejTotal = 0
    nIter = 0
    while nIter < nMaxIter:
        stdDev = np.std(data)
        mean = np.mean(data)
        oldSize = data.size
        data = data[(data <= mean + nSigma * stdDev) &
                    (data >= mean - nSigma * stdDev)]
        nSrej = oldSize - data.size
        nSrejTotal += nSrej
        if nSrej == 0:
            break
        else:
            nIter += 1
    return data, nSrejTotal


def calcCentroid(data, width=3., binSize=0.1):
    """
    Calculates the intensity averaged centroid of a set of data
    """

    mean = np.mean(data)
    stdDev = np.std(data)

    if stdDev == 0.:
        return mean, 0.0

    width *= stdDev
    binSize *= stdDev

    bins = np.arange(mean-width, mean+width+binSize, binSize)

    hist = np.histogram(data, bins=bins, range=(mean-width, mean+width))
    midPoints = np.array([np.mean([hist[1][ii], hist[1][ii+1]]) for ii in range(0, hist[1].size-1)])

    return np.sum(hist[0] * midPoints) / np.sum(hist[0]), stdDev


def calcSky(data, xCentre, yCentre, annulus, dannulus, refine=False, nRefine=10):
    """
    Calculates the sky value around a certain position given an annulus (radii)
    and dannulus (width of the annulus).
    """

    dataReg, mask = calcMask(data, xCentre, yCentre, annulus,
                             r2=annulus+dannulus, refine=False)
    dataAnn = dataReg * mask
    dataAnn = dataAnn[np.isfinite(dataAnn)]

    validData, nSrej = getRejPixels(dataAnn)
    nSky = validData.size

    mSky, stdDev = calcCentroid(validData)

    del dataReg, dataAnn, validData

    # Returns value of sky (per pixel), standard deviation, number of pixels
    # used for the sky estimation and number of pixels rejected
    return mSky, stdDev, nSky, nSrej


def phot(data, centres, apertures, annulus=None, dannulus=None, logFile=None,
         ext=0, refine=False, nRefine=10):
    """
    Calculates the aperture photometry with sky subtraction around centres for
    a number of apertures.

    Output:
        returnFluxes = [[centreID, xCentre, yCentre, aperture, sum, area, flux], ...]
        returnSky = [[mSky, stdDev, nSky, nSrej], ... (for each centre)]
    """

    # If data is a string, read the file, otherwise, data is supposed to be
    # a numpy ndarray
    if isinstance(data, str):
        import astropy.io.fits as pf
        data = pf.getdata(data, ext)

    if not isinstance(data, np.ndarray):
        raise TypeError('The required argument \'data\' is not a numpy.ndarray.')

    centres = np.array(centres)
    apertures = np.array(apertures)
    nCentres = centres.shape[0]
    nApertures = apertures.size
    returnFluxes = np.zeros((nApertures * nCentres, 7), dtype=np.float)
    returnSky = np.zeros((nCentres, 4), dtype=np.float)

    # If logFile, opens a unit for output
    if logFile is not None:
        unitOut = open(logFile, 'w')

    for nn in range(centres.shape[0]):

        xCentre = centres[nn, 0]
        yCentre = centres[nn, 1]

        # If annulus and dannulus are defined, calculates the
        # sky background
        if annulus is not None and dannulus is not None:
            try:
                annulus = float(annulus)
                dannulus = float(dannulus)
            except:
                raise ValueError('Annulus and/or dannulus are not numeric.')

            mSky, stdDev, nSky, nSrej = calcSky(data, xCentre, yCentre, annulus, dannulus,
                                                refine=refine, nRefine=nRefine)

        else:

            mSky = stdDev = nSky = nSrej = np.nan

        sums = np.zeros(nApertures)
        fluxes = np.zeros(nApertures)
        areas = np.zeros(nApertures)

        for ii in range(nApertures):

            dataReg, mask = calcMask(data, xCentre, yCentre, apertures[ii], refine=refine, nRefine=nRefine)
            dataAp = dataReg * mask
            dataAp = dataAp[np.isfinite(dataAp)]

            sums[ii] = np.sum(dataAp)

            if refine is None:
                areas[ii] = dataAp.size
            else:
                areas[ii] = np.pi * apertures[ii]**2

            if mSky is not np.nan:
                fluxes[ii] = sums[ii] - areas[ii] * mSky
            else:
                fluxes[ii] = sums[ii]

        returnFluxes[nn*nApertures:(nn+1)*nApertures, :] = np.array([[nn+1]*nApertures,
                                                                     [xCentre]*nApertures,
                                                                     [yCentre]*nApertures,
                                                                     apertures, sums,
                                                                     areas, fluxes]).T
        returnSky[nn, :] = (mSky, stdDev, nSky, nSrej)

        del dataReg, dataAp

        # Outputs the data
        if logFile is not None:

            print >>unitOut, '# %12s %12s' % ('XCENTRE', 'YCENTRE')
            print >>unitOut, '# %12.3E %12.3E' % (xCentre, yCentre)

            if annulus is not None and dannulus is not None:
                print >>unitOut, '# %12s %12s %12s %12s' % ('MSKY', 'STDEV', 'NSKY', 'NSREJ')
                print >>unitOut, '# %12.3E %12.3E %12.3E %12.3E' % (mSky, stdDev, nSky, nSrej)
            else:
                print >>unitOut, '# %12s %12s %12s %12s' % ('MSKY', 'STDEV', 'NSKY', 'NSREJ')
                print >>unitOut, '# %12s %12s %12s %12s' % ('N/A', 'N/A', 'N/A', 'N/A')

            writer = at.FixedWidth if nn == 0 else at.FixedWidthNoHeader

            at.write([nApertures*[xCentre], nApertures*[yCentre],
                     apertures, sums, areas, fluxes], unitOut,
                     Writer=writer, names=['XCENTRE', 'YCENTRE', 'AP', 'SUM', 'AREA', 'FLUX'])

            print >>unitOut

    del data

    return returnFluxes, returnSky


def photTest():
    """
    Testing routine
    """

    print

    image = './TestImages/StdImage.fits'
    centres = np.loadtxt('./TestImages/coordsStd.coo')
    # centres = np.array([[501, 501]])
    apertures = np.arange(1, 16)

    phot(image, centres, apertures,
         annulus=15, dannulus=2, refine=True, nRefine=10,
         logFile='test.log')

    print
    return


# if __name__ == '__main__':
#     photTest()
