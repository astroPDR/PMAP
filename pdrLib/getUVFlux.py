#!/usr/bin/env python
# encoding: utf-8
"""
getUVFlux.py

Created by José Sánchez-Gallego on 3 May 2013.
Copyright (c) 2013. All rights reserved.

FUV flux measurements
First, try to read the flux file. Otherwise, calculate the fluxes. Output: fluxtable

"""

import os
import numpy as np


def getUVFlux(configOpts, logger):

    print "Measuring FUV fluxes ..."
    uvheader, uvimage = pdr.open_image(configOpts['fuvImage'])
    data_uv = configOpts['data_uv']
    coords_uv = os.path.splitext(configOpts['fuvImage'])[0] + '_Peaks.dat' #as defined earlier

    try:
        dummy = file(data_uv, 'r')
        print "... fluxes read from file {0} instead.".format(data_uv)
        columns = np.dtype([('PDRID', 'int',), ('RA', '|S14'), ('DEC', '|S14'), ('aperture', 'float'), ('mean_at_r', 'float'), ('bgflux', 'float'), ('cumulflux', 'float'), ('netflux', 'float'), ('sflux', 'float')])
        fluxtable = pdr.read_array(data_uv, columns)
        print "  {0} records read".format(np.size(fluxtable.PDRID))
    except:
        #file does not exist, so we calculate. Force recalculation by removing the file.
        #Variables in: uvheader, uvimage, uvcoords
        #Variable out: fluxtable (includes uvcoords afterwards)
        #   fluxtable labels: PDRID, RA, DEC, aperture, mean_at_r, bgflux, cumulflux, netflux, sflux
        print "No previous flux measurements."
        uvcoords = pdr.read_coords(coords_uv)
        print "Flux will be measured at {0} positions on image {1} ...".format(len(uvcoords), configOpts['fuvImage'])
        print "Fluxes will be scaled by a factor {0}".format(configOpts['FUVscale'])
        #This scaling is useful since clumpfind doesn't seem to work properly with ~1e-15 values.
        fluxtable = pdr.fuv_flux(configOpts['fuvImage'], uvheader, uvcoords, configOpts, verbose = True)
        print "... done"

        #print fluxtable.RA, fluxtable.RA[0] #here fluxtable.RA[0] is simply a string

        #Write: UV fluxes
        fluxfile = fancyPrint(writeLog=True, logFile=data_uv, verbose=True)
        print "Flux info:"
        #print "#PDRID, aperture (arcsec), mean at r (units/arcsec^2), bg flux, cumul. flux, net flux"
        fluxfile.write("#PDRID, RA, DEC, aperture (arcsec), mean at r (units/pixel), bg flux, cumul. flux, net flux, sigma net flux\n")
        for n in range(len(fluxtable)):
            #fedora bug: changed PDRID:>3 to PDRID:3g
            #print "{0.PDRID:3g}, {1.aperture:5.2f}, {1.mean_at_r:7.5e}, {1.bgflux:7.5e}, {1.cumulflux:7.5e}, {1.netflux:7.5e}".format(uvcoords[n], fluxtable[n])
            fluxfile.write("{0.PDRID:3g}, {0.RA}, {0.DEC}, {1.aperture:7.5e}, {1.mean_at_r:7.5e}, {1.bgflux:7.5e}, {1.cumulflux:7.5e}, {1.netflux:7.5e}, {1.sflux:7.5e}\n".format(uvcoords[n], fluxtable[n]))
          #fluxfile.close()

    print
    #print "Exiting here."
    #exit(0)

    return uvheader, fluxtable
