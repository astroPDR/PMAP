#!/usr/bin/env python
# encoding: utf-8
'''
getUVFlux.py

Created by José Sánchez-Gallego on 3 May 2013.
Copyright (c) 2013. All rights reserved.

FUV flux measurements
First, try to read the flux file. Otherwise, calculate the fluxes. Output: fluxtable

'''

import os
from .. import at, pf, pw
from fluxes import fuv_flux


def getUVFlux(configOpts, logger):

    logger.write('Measuring FUV fluxes ...', newLine=True)
    data_uv = configOpts['data_uv']

    if os.path.exists(data_uv) and configOpts['overwrite'] is False:
        try:
            # columns = np.dtype(
            #     [('PDRID', 'int',), ('RA', '|S14'), ('DEC', '|S14'),
            #      ('aperture', 'float'), ('mean_at_r', 'float'),
            #      ('bgflux', 'float'), ('cumulflux', 'float'),
            #      ('netflux', 'float'), ('sflux', 'float')])
            # fluxtable = read_array(data_uv, columns)
            fluxTable = at.read(data_uv, delimiter=' ')
            logger.write(
                '... fluxes read from file {0} instead.'.format(data_uv))
            logger.write(
                '... {0} records read'.format(len(fluxTable)))

            return fluxTable

        except:
            pass

    # file does not exist, so we calculate. Force recalculation by removing the file.
    # Variables in: uvheader, uvimage, uvcoords
    # Variable out: fluxtable (includes uvcoords afterwards)
    # fluxtable labels: PDRID, RA, DEC, aperture, mean_at_r, bgflux,
    # cumulflux, netflux, sflux

    fuvImage = configOpts['fuvImage']
    hduFUV = pf.open(fuvImage)
    fuvWCS = pw.WCS(hduFUV[0].header)
    fuvPeaksFile = configOpts['fuvPeaksFile']

    logger.write('No previous flux measurements.')
    fuvCoords = at.read(fuvPeaksFile, delimiter=',')

    logger.write(
        'Flux will be measured at {0} positions on image {1} ...'.format(len(fuvCoords),
                                                                         fuvImage))
    logger.write(
        'Fluxes will be scaled by a factor {0}'.format(configOpts['FUVscale']))
    # This scaling is useful since clumpfind doesn't seem to work properly
    # with ~1e-15 values.
    fluxTable = fuv_flux(fuvImage, fuvWCS, fuvCoords, configOpts, logger)

    # logger.write(fluxtable.RA, fluxtable.RA[0] #here fluxtable.RA[0] is
    # simply a string

    # Write: UV fluxes
    at.write(fluxTable, data_uv, format='fixed_width', delimiter=' ')

    # fluxfile = open(data_uv, 'w')
    # # logger.write('#PDRID, aperture (arcsec), mean at r (units/arcsec^2),
    # # bg flux, cumul. flux, net flux'
    # fluxfile.write(
    #     '#PDRID, RA, DEC, aperture (arcsec), mean at r (units/pixel), ' +
    #     'bg flux, cumul. flux, net flux, sigma net flux\n')
    # for n in range(len(fluxTable)):
    #     # fedora bug: changed PDRID:>3 to PDRID:3g
    #     # logger.write('{0.PDRID:3g}, {1.aperture:5.2f}, {1.mean_at_r:7.5e}, ' +
    #     #              '{1.bgflux:7.5e}, {1.cumulflux:7.5e}, ' +
    #     #              '{1.netflux:7.5e}'.format(uvcoords[n], fluxtable[n])
    #     fluxfile.write('{0.PDRID:3g}, {0.RA}, {0.DEC}, {1.aperture:7.5e}, ' +
    #                    '{1.mean_at_r:7.5e}, {1.bgflux:7.5e}, ' +
    #                    '{1.cumulflux:7.5e}, {1.netflux:7.5e}, ' +
    #                    '{1.sflux:7.5e}\n'.format(uvcoords[n], fluxTable[n]))
    # fluxfile.close()

    return fluxTable
