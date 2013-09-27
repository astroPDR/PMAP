#!/usr/bin/env python
# encoding: utf-8
"""
fluxes.py

Created by José Sánchez-Gallego on 8 May 2013.
Copyright (c) 2013. All rights reserved.

This file contains a series of functions intended to measure fluxes.

"""


def generic_flux(filename, coordsPix, maxrad_pix, pix_size = 1.0, contrast = 0.5, fluxerror = 0.3, logfile = False):

    from calcFlux import calcFluxMulti

    genericFluxResults = [] # The list with the final results to return

    # Calls calcFluxMulti an gets the areas and fluxes for each position and aperture
    # from 1 to maxrad_pix
    fluxes = calcFluxMulti(filename, coordsPix, 1, maxrad_pix, step=1)

    for pos in range(coordsPix.shape[0]): # iterates over the different positions

        fluxlist = []
        avglist = []
        mean_at_r = 0; aperture = 0; bgflux = 0; cumulflux = 0; netflux = 0; sflux = 0

        for n in range(maxrad_pix): #iterate over rings of width 1 pixel

            r = n + 1
            fluxlist.append([fluxes[pos][5][n], fluxes[pos][6][n]])
            # print fluxlist
            if n > 0:
                #incremental flux divided by ring surface area
                ring_avg = (fluxlist[n][1] - fluxlist[n-1][1]) / (fluxlist[n][0] - fluxlist[n-1][0])
                avglist.append(ring_avg) #current avglist index will be n-1 since it starts at n=1

                if (logfile):
                    logdump(logfile, "{0:3g} {1:7.2e} {2:8.3e} {3:7.5e}\n".format(r, fluxlist[n][0], fluxlist[n][1], ring_avg))

            if (n > 1) and (not mean_at_r):
                #Check if the mean value has dropped to at least 50% (contrast) of the initial value (arbitrary)
                #  note that if the actual max value is much higher, 50% is not indicative of any kind of source contrast
                #Then check if a local minimum was encountered.
                #This could create problems with very poor source contrast, but those may need to be filtered out earlier on
                #For example: need to drop to at least half the initial mean value first
                if (avglist[n-1] < (avglist[0] * contrast)):
                    if (avglist[n-1] > avglist[n-2]): #list is 1 shorter than n
                        mean_at_r = avglist[n-2] #per pixel
                        aperture = (r-1) * pix_size #in arcsec and assuming square pixels
                        bgflux = mean_at_r * fluxlist[n-2][0] #flux/pixel*area[in pixels]
                        cumulflux = fluxlist[n-2][1]
                        netflux = cumulflux - bgflux
                        mean_at_r = mean_at_r / (pix_size*pix_size) #!need to report per arcsec^2
                        sflux = fluxerror*netflux #flux error as in config file (fixed relative error)

                        if (logfile):
                          logdump(logfile, "#First local minimum registered\n")
                        #print aperture, r-1

        if (not mean_at_r): #we exit leaving all the values at zero
            print "Warning: minimum not found! (pos {0})".format(pos)
            if (logfile):
                logdump(logfile, "#Warning: minimum not found!\n")

        genericFluxResults.append([mean_at_r, aperture, bgflux, cumulflux, netflux, sflux])

    return genericFluxResults


# fuv_flux
def fuv_flux(imageName, header, coords, configOpts, verbose=False):
  #Wishlist: use qphot to get all the annuli in one go
  #configOpts used: FUVscale, contrast, maxrad, uv_log
  #Determines the (FUV) flux at a specified coordinate
  #A net flux is calculated based on the a local minimum in ring-averaged flux values per pixel
  #The maximum radius can't be bigger than half the box size (not checked)
  #Need to implement: option to accept first local minimum or go to next

  logfile = configOpts['uv_log']

  columns = np.dtype([('PDRID', 'int32'), ('RA', '|S14'), ('DEC', '|S14'), ('aperture', 'float'), ('mean_at_r', 'float'), ('bgflux', 'float'), ('cumulflux', 'float'), ('netflux', 'float'), ('sflux', 'float')])
  results = []#[[] for dummy in xrange(len(columns))]
  pix_size = [header.getXPixelSizeDeg() * 3600, header.getYPixelSizeDeg() * 3600]
  #pix_area = pix_size[0] * pix_size[1]
  logdump(logfile, "#Number of coordinates: {0}\n".format(len(coords)))
  logdump(logfile, "#R(pix) Area (pix) Flux (image units) Flux/pix\n")

  #Determine flux for all coordinates

  maxrad_pix = int(round(configOpts['maxrad']/pix_size[0])) # Global maximum radius (in pixels)
  coordsPix = np.zeros([len(coords), 2], np.float)

  for p in range(len(coords)):
    #Define circle center and check pixel size
    x, y = np.round(header.wcs2pix(astCoords.hms2decimal(strip(coords[p].RA), " "), astCoords.dms2decimal(strip(coords[p].DEC), " "))) #x, y not integers so int(x) will trunc not round
    #x, y should be rounded values or requesting dummy[y,x] will lead to truncation errors
    logdump(logfile, "#Source {0} at pixel position: {1}, {2}\n".format(coords[p].PDRID, x, y))
    coordsPix[p, :] = [x, y]

  genericFluxResults =  generic_flux(imageName, coordsPix, maxrad_pix, pix_size[0],
    configOpts['contrast'], configOpts['fluxerror'], logfile = logfile)

  for p in range(len(coords)):

    mean_at_r, aperture, bgflux, cumulflux, netflux, sflux = genericFluxResults[p]

    #Scale the results -- for generic use set FUVscale to 0
    scale = configOpts['FUVscale'] # Converts from cps to flux units
    mean_at_r *= scale
    bgflux *= scale
    cumulflux *= scale
    netflux *= scale
    sflux *= scale

    #Store the results
    if verbose:
      print "local minimum at", mean_at_r, aperture
      print "bgflux", bgflux
      print "cumul flux", cumulflux
      print "net flux", netflux

    results.append([coords[p].PDRID, coords[p].RA, coords[p].DEC, aperture, mean_at_r, bgflux, cumulflux, netflux, sflux])

  return np.rec.array(results, dtype=columns)

