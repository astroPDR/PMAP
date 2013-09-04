#!/usr/bin/env python
# encoding: utf-8
"""
phot.py

Created by José Sánchez-Gallego on 3 May 2013.
Copyright (c) 2013. All rights reserved.

"""

import numpy as np
from scipy.stats.stats import nanmean
from .. import table
from .. import at


def rebin(aa, shape):
    """
    Resizes an array to a certain shape
    """
    sh = shape[0], aa.shape[0] // shape[0], shape[1], aa.shape[1] // shape[1]
    return nanmean(nanmean(aa.reshape(sh), axis=-1), axis=1)


def calcMask(xCentre, yCentre, shape, r1, r2=None, refine=False, nRefine=10):
    """
    Returns a mask with the pixels included in the aperture, as well as a list of the
    pixels in the original image. If r2 is defined, the apertures is annular.
    If refine, the data is re-meshed by a factor nRefine.
    """

    edge = 0

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
    ii = yCentre + np.arange(-rMax, rMax+1, step)
    jj = xCentre + np.arange(-rMax, rMax+1, step)

    # Takes care of edges of the image
    ii[(ii < 0) | (ii > shape[0]-1)] = np.nan
    jj[(jj < 0) | (jj > shape[1]-1)] = np.nan

    # If it touches the edge of the image and this is not an annular aperture, fails.
    if True in np.isnan(ii) | np.isnan(jj):
        edge = 1
        if r2 is None:
            return None, None, edge

    # Coords grid
    coords = np.dstack(np.meshgrid(ii, jj)).reshape(-1, 2)

    # Calculates distances
    distances = ((coords[:, 0]-yCentre)**2 + (coords[:, 1]-xCentre)**2).reshape(ii.shape[0], -1)
    mask = np.zeros(distances.shape, np.float)

    # Determines pixels within aperture and set them to 1.
    if r2 is None:
        mask[distances <= (r1**2)] = 1
    else:
        mask[(distances >= (r1**2)) & (distances <= (r2**2))] = 1

    # Rebins the mask to its original size
    if step < 1.:
        mask = rebin(mask, (ii.shape[0] // nRefine, jj.shape[0] // nRefine))
        ii = yCentre + np.arange(-rMax, rMax+1)
        jj = xCentre + np.arange(-rMax, rMax+1)
        coords = np.dstack(np.meshgrid(ii, jj)).reshape(-1, 2)

    # To avoid memory usagem we return only those pixels within the aperture
    mask = mask.flatten()
    nonZeroMask = np.where(mask != 0.)
    mask = mask[nonZeroMask]
    coords = coords[nonZeroMask]

    return (coords[:, 0].astype(int), coords[:, 1].astype(int)), mask, edge


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

    coords, mask, edge = calcMask(xCentre, yCentre, data.shape, annulus,
                                  r2=annulus+dannulus, refine=False)

    dataAnn = data[coords] * mask

    validData, nSrej = getRejPixels(dataAnn)
    nSky = validData.size

    mSky, stdDev = calcCentroid(validData)

    del dataAnn, validData

    # Returns value of sky (per pixel), standard deviation, number of pixels
    # used for the sky estimation and number of pixels rejected
    return mSky, stdDev, nSky, nSrej, edge


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
    # nCentres = centres.shape[0]
    nApertures = apertures.size
    returnData = []
    # returnFluxes = np.zeros((nApertures * nCentres, 8), dtype=np.float)
    # returnSky = np.zeros((nCentres, 5), dtype=np.float)

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

            mSky, stdDev, nSky, nSrej, edgeAnn = calcSky(data, xCentre, yCentre, annulus, dannulus,
                                                         refine=refine, nRefine=nRefine)

        else:

            mSky = stdDev = nSky = nSrej, edgeAnn = np.nan

        sums = np.zeros(nApertures)
        fluxes = np.zeros(nApertures)
        areas = np.zeros(nApertures)
        edges = np.zeros(nApertures, np.int)

        for ii in range(nApertures):

            coords, mask, edgeAp = calcMask(xCentre, yCentre, data.shape, apertures[ii],
                                            refine=refine, nRefine=nRefine)

            if edgeAp is 0:

                dataAp = data[coords] * mask

                sums[ii] = np.sum(dataAp)

                if refine is None:
                    areas[ii] = dataAp.size
                else:
                    areas[ii] = np.pi * apertures[ii]**2

                if mSky is not np.nan:
                    fluxes[ii] = sums[ii] - areas[ii] * mSky
                else:
                    fluxes[ii] = sums[ii]

                del dataAp

            else:

                sums[ii] = areas[ii] = fluxes[ii] = np.nan
                edges[ii] = 1

        # returnFluxes[nn*nApertures:(nn+1)*nApertures, :] = np.array([[nn+1]*nApertures,
        #                                                              [xCentre]*nApertures,
        #                                                              [yCentre]*nApertures,
        #                                                              apertures, sums,
        #                                                              areas, fluxes, edges]).T
        # returnSky[nn, :] = (mSky, stdDev, nSky, nSrej, edgeAnn)

        # Outputs the data
        if logFile is not None:

            print >>unitOut, '# %12s %12s' % ('XCENTRE', 'YCENTRE')
            print >>unitOut, '# %12.3E %12.3E' % (xCentre, yCentre)

            if annulus is not None and dannulus is not None:
                print >>unitOut, '# %12s %12s %12s %12s %12s' % \
                    ('MSKY', 'STDEV', 'NSKY', 'NSREJ', 'EDGE')
                print >>unitOut, '# %12.3E %12.3E %12.3E %12.3E %12d' % \
                    (mSky, stdDev, nSky, nSrej, edgeAnn)
            else:
                print >>unitOut, '# %12s %12s %12s %12s %12s' % \
                    ('MSKY', 'STDEV', 'NSKY', 'NSREJ', 'EDGE')
                print >>unitOut, '# %12s %12s %12s %12s %12s' % \
                    ('N/A', 'N/A', 'N/A', 'N/A', 'N/A')

            writer = at.FixedWidth if nn == 0 else at.FixedWidthNoHeader

            at.write(
                [nApertures*[xCentre], nApertures*[yCentre],
                 apertures, sums, areas, fluxes, edges],
                unitOut, Writer=writer,
                names=['XCENTRE', 'YCENTRE', 'AP', 'SUM', 'AREA', 'FLUX', 'EDGE'])

            print >>unitOut

        # Creates the Table that will be included in the output list
        edges = np.array(edges, np.bool)
        tmpTable = table.Table([apertures, sums, areas, fluxes, edges],
                               names=['aperture', 'sum', 'area', 'flux', 'edge'],
                               meta={'id': nn+1, 'xCentre': xCentre, 'yCentre': yCentre,
                                     'mSky': mSky, 'stdDev': stdDev, 'nSky': nSky,
                                     'nSrej': nSrej, 'edgeSky': bool(edgeAnn)})
        returnData.append(tmpTable)
        del tmpTable

    del data
    if logFile is not None:
        unitOut.close()

    return returnData


# def photTest():
#     """
#     Testing routine
#     """

#     print

#     image = './TestImages/StdImage.fits'
#     centres = np.loadtxt('./TestImages/coordsStd.coo')
#     # centres = np.array([[501, 501]])
#     apertures = np.arange(1, 16)

#     phot(image, centres, apertures,
#          annulus=15, dannulus=2, refine=True, nRefine=10,
#          logFile='test.log')

#     print
#     return


# if __name__ == '__main__':
#     photTest()
