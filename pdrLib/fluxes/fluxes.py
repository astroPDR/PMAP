#!/usr/bin/env python
# encoding: utf-8
"""
fluxes.py

Created by José Sánchez-Gallego on 8 May 2013.
Copyright (c) 2013. All rights reserved.

This file contains a series of functions intended to measure fluxes.

"""

import numpy as np
from phot import phot
from .. import Error
from .. import table


def generic_flux(filename, coordsPix, maxrad_pix, logger,
                 pix_size=1.0, contrast=0.5, fluxerror=0.3,
                 refine=True, nRefine=10):

    genericFluxResults = []  # The list with the final results to return

    # Calls calcFluxMulti an gets the areas and fluxes for each position and aperture
    # from 1 to maxrad_pix
    logger.write('Measuring fluxes now ...', doLog=False)
    fluxes = phot(filename, coordsPix, np.arange(1, maxrad_pix+1),
                  annulus=maxrad_pix, dannulus=1, refine=refine, nRefine=nRefine,
                  logFile=None)

    for pos in range(coordsPix.shape[0]):  # iterates over the different positions

        fluxlist = []
        avglist = []
        mean_at_r = 0
        aperture = 0
        bgflux = 0
        cumulflux = 0
        netflux = 0
        sflux = 0

        fluxPos = fluxes[pos]
        if fluxPos.meta['id'] != pos+1:
            Error.raiseError('Error matching phot measurement with region')

        for n in range(maxrad_pix):  # iterate over rings of width 1 pixel

            r = n + 1
            fluxlist.append([fluxPos['area'][n], fluxPos['flux'][n]])
            # print fluxlist
            if n > 0:
                # incremental flux divided by ring surface area
                ring_avg = (fluxlist[n][1] - fluxlist[n-1][1]) / (fluxlist[n][0] -
                                                                  fluxlist[n-1][0])
                # current avglist index will be n-1 since it starts at n=1
                avglist.append(ring_avg)

                # logger.write('{0:3g} {1:7.2e} {2:8.3e} {3:7.5e}'.format(r, fluxlist[n][0],
                #                                                         fluxlist[n][1],
                #                                                         ring_avg),
                #              doPrint=True)

            if (n > 1) and (not mean_at_r):
                """
                Check if the mean value has dropped to at least 50% (contrast) of
                   the initial value (arbitrary) note that if the actual max value is
                   much higher, 50% is not indicative of any kind of source contrast
                Then check if a local minimum was encountered.
                This could create problems with very poor source contrast, but those may
                   need to be filtered out earlier on
                For example: need to drop to at least half the initial mean value first
                """

                if (avglist[n-1] < (avglist[0] * contrast)):
                    if (avglist[n-1] > avglist[n-2]):          # list is 1 shorter than n
                        mean_at_r = avglist[n-2]               # per pixel
                        aperture = (r-1) * pix_size            # in arcsec and assuming square pix
                        bgflux = mean_at_r * fluxlist[n-2][0]  # flux/pixel*area[in pixels]
                        cumulflux = fluxlist[n-2][1]
                        netflux = cumulflux - bgflux
                        mean_at_r = mean_at_r / (pix_size*pix_size)  # need to report per arcsec^2
                        # flux error as in config file (fixed relative error)
                        sflux = fluxerror*netflux

                        # logger.write('First local minimum registered', doPrint=True)
                        # print aperture, r-1

        if (not mean_at_r):  # we exit leaving all the values at zero
            logger.write('Warning: minimum not found! (pos {0})'.format(pos))
            # logger.write('Warning: minimum not found!', doPrint=True)

        genericFluxResults.append([mean_at_r, aperture, bgflux, cumulflux, netflux, sflux])

    return genericFluxResults


# fuv_flux
def fuv_flux(imageName, header, coords, configOpts, logger):
    """
    configOpts used: FUVscale, contrast, maxrad, uv_log
    Determines the (FUV) flux at a specified coordinate
    A net flux is calculated based on the a local minimum
       in ring-averaged flux values per pixel
    The maximum radius can't be bigger than half the box size (not checked)
    Need to implement: option to accept first local minimum or go to next
    """

    results = []  # [[] for dummy in xrange(len(columns))]
    #pix_area = pix_size[0] * pix_size[1]

    pix_size = configOpts['FUV_Pix_Scale']

    # logger.write('Number of coordinates: {0}'.format(len(coords)), newLine=True, doPrint=True)

    # Determine flux for all coordinates
    maxrad_pix = int(round(configOpts['maxrad'] / pix_size))  # Global maximum radius (in pix)
    coordsPix = np.zeros([len(coords), 2], np.float)

    for p in range(len(coords)):
        # Define circle center and check pixel size
        # x, y not integers so int(x) will trunc not round
        # worldCoord = np.array([[astCoords.hms2decimal(strip(coords[p]['RA']), ' '),
                               # astCoords.dms2decimal(strip(coords[p]['DEC']), ' ')]])
        # x, y = np.round(header.wcs_world2pix(worldCoord, 0)[0])
        x = int(np.round(coords[p]['XPIXEL']-1))
        y = int(np.round(coords[p]['YPIXEL']-1))

        # x, y should be rounded values or requesting dummy[y,x] will lead to truncation errors
        # logger.write('Source {0} at pixel position: {1}, {2}'.format(coords[p]['PDRID'], x, y),
                     # doPrint=False)
        coordsPix[p, :] = [x, y]

    # logger.write('R (pix) Area (pix) Flux (image units) Flux/pix', doPrint=True)

    genericFluxResults = generic_flux(imageName, coordsPix, maxrad_pix, logger,
                                      pix_size=pix_size, contrast=configOpts['contrast'],
                                      fluxerror=configOpts['fluxerror'],
                                      refine=configOpts['refineFlux'],
                                      nRefine=configOpts['refineFactor'])

    for p in range(len(coords)):

        mean_at_r, aperture, bgflux, cumulflux, netflux, sflux = genericFluxResults[p]

        # Scale the results -- for generic use set FUVscale to 0
        scale = configOpts['FUVscale']  # Converts from cps to flux units
        mean_at_r *= scale
        bgflux *= scale
        cumulflux *= scale
        netflux *= scale
        sflux *= scale

        # Store the results
        # if configOpts['verbose'] is True:
        #     logger.write('local minimum at', mean_at_r, aperture, doPrint=False, doLog=True)
        #     logger.write('bgflux', bgflux, doPrint=False, doLog=True)
        #     logger.write('cumul flux', cumulflux, doPrint=False, doLog=True)
        #     logger.write('net flux', netflux, doPrint=False, doLog=True)

        results.append([coords[p]['PDRID'], coords[p]['RA'], coords[p]['DEC'],
                        aperture, mean_at_r, bgflux, cumulflux, netflux, sflux])

    names = ('PDRID', 'RA', 'Dec', 'aperture', 'mean_at_r', 'bgflux',
             'cumulflux', 'netflux', 'sflux')

    return table.Table(np.rec.array(results, names=names))
