#!/usr/bin/env python
# encoding: utf-8
"""
createCatalog.py

Created by José Ramón Sánchez-Gallego on 2011-03-25.
JSH 2011-10-17: updated to deal with negative DEC correctly

JRSG 2013-04-22.
    -- General clean-up.
    -- Now uses logger.
    -- Dependences are now managed using impPyFits, impPyWCS and impVO.
    -- If vo is not available, a FITS table is created.

Copyright (c) 2011. All rights reserved.

This script creates a VOTable catalog with the properties of each region within
the selected mask.

"""

import numpy as np
import os
from fitellipse import fitellipse
from Regions import *
import colorama as cm

from impPyFits import *
from impPyWCS import *
from impVO import *

from Error import raiseWarning, raiseError


# Main routine which read the image and mask files, calculates various parameters and
# creates the VOTable catalogue and the DS9 regions files.
def createCatalog(image, mask, voTableCat, ds9RegsFile, logger,
                  peaksFile=None, ellipse=False, ellMode='linear', plot=False):

    if voModule is False:
        logger.write(cm.Fore.RED + 'No VO module found. A FITS table will be generated instead.' +
                     cm.Style.RESET_ALL, newLine=True)

# Reads the header, data and loads the astrometric solution
    hduImage = pf.open(image)
    hduMask = pf.open(mask)
    try:
        wcsMask = pw.WCS(hduMask[0].header)
    except:
        wcsMask = None
        logger.write(cm.Fore.YELLOW + 'No WCS information. The VOTable file will not include those fields.' +
                     cm.Style.RESET_ALL)

    logger.write('Gathering mask data ... ', newLine=True)
    regions = RegionSet(hduMask[0].data, image=hduImage[0].data)

    logger.write('Preparing data ... ')
    centroids = np.array([regions.Regions[idReg].getCentroid()[0] for idReg in regions.idRegions])
    radii = [regions.Regions[idReg].getCentroid()[1] for idReg in regions.idRegions]
    sizes = [regions.Regions[idReg].where[0].size for idReg in regions.idRegions]
    peaks = [regions.Regions[idReg].peak for idReg in regions.idRegions]
    peaksPix = np.array([regions.Regions[idReg].peakWhere for idReg in regions.idRegions]) + 1
    fluxes = [regions.Regions[idReg].getFlux() for idReg in regions.idRegions]

    if wcsMask is not None:
        skyCoordsCentroids = wcsMask.wcs_pix2sky(np.fliplr(centroids), 0)
        skyCoordsPeak = wcsMask.wcs_pix2sky(np.fliplr(peaksPix), 0)

        try:
            scale = np.abs(hduImage[0].header['CD1_1']) * 3600.
        except:
            try:
                scale = np.abs(hduImage[0].header['CDELT1']) * 3600.
            except:
                wcsMask = None
                logger.write(cm.Fore.YELLOW + 'No WCS information. The VOTable file will ' +
                             'not include those fields.' + cm.Style.RESET_ALL)

        radiiSky = radii * scale

# Calculates an ellipse fit for the pixels in the region
    if ellipse:
        logger.write('Fitting ellipse parameters ... ')

        def getEllParams(data, id):
            try:
                return list(fitellipse(data, ellMode)) + [0]
            except RuntimeError:
                return [regions.Regions[id].getCentroid()[0][0], regions.Regions[id].getCentroid()[0][1],
                        regions.Regions[idReg].getCentroid()[1], regions.Regions[idReg].getCentroid()[1],
                        0.0, 1]

        ellipseFit = np.array([getEllParams(regions.Regions[idReg].pixels, idReg) for idReg in regions.idRegions])
        ellipseCentroidWCS = np.array(wcsMask.wcs_pix2sky(ellipseFit[:, [1, 0]], 0))
        ellipseAxisWCS = np.array(ellipseFit[:, [2, 3]] * scale)
        ellipseAngle = np.array(-ellipseFit[:, 4] * 180. / np.pi + 90.)

    if voModule is True:
        # Creates a VOTable with the data
        logger.write('Creating VOTable ... ')

        if os.path.exists(voTableCat):
            os.remove(voTableCat)

        # from vo.tree import VOTableFile, Resource, Table, Field

        votable = VOTableFile()
        resource = Resource()
        votable.resources.append(resource)
        table = Table(votable)
        resource.tables.append(table)

        # Define fields
        fieldsVO = [Field(votable, name='ID', datatype='int'),
                    Field(votable, name='X_CENTROID', datatype='float', unit='pix'),
                    Field(votable, name='Y_CENTROID', datatype='float', unit='pix'),
                    Field(votable, name='RADIUS_PIX', datatype='float', unit='pix'),
                    Field(votable, name='SIZE', datatype='float', unit='pix')]

        if wcsMask is not None:
            fieldsVO += [Field(votable, name='RA_CENTROID', datatype='double', unit='deg'),
                         Field(votable, name='DEC_CENTROID', datatype='double', unit='deg'),
                         Field(votable, name='RADIUS_SKY', datatype='float', unit='arcsec')]

        fieldsVO += [Field(votable, name='PEAK', datatype='float', unit='CPS'),
                     Field(votable, name='X_PEAK', datatype='float', unit='pix'),
                     Field(votable, name='Y_PEAK', datatype='float', unit='pix')]

        if wcsMask is not None:
            fieldsVO += [Field(votable, name='RA_PEAK', datatype='double', unit='deg'),
                         Field(votable, name='DEC_PEAK', datatype='double', unit='deg')]

        fieldsVO += [Field(votable, name='FLUX', datatype='float', unit='CPS')]

        if ellipse is True:
            fieldsVO += [Field(votable, name='X_ELLIPSE', datatype='float', unit='pix'),
                         Field(votable, name='Y_ELLIPSE', datatype='float', unit='pix'),
                         Field(votable, name='SEMIAX1_PIX', datatype='float', unit='pix'),
                         Field(votable, name='SEMIAX2_PIX', datatype='float', unit='pix'),
                         Field(votable, name='ELL_ANGLE', datatype='float', unit='deg'),
                         Field(votable, name='ELL_ERROR', datatype='double')]

            if wcsMask is not None:
                fieldsVO += [Field(votable, name='RA_ELLIPSE', datatype='double', unit='deg'),
                             Field(votable, name='DEC_ELLIPSE', datatype='double', unit='deg'),
                             Field(votable, name='SEMIAX1_SKY', datatype='float', unit='arcsec'),
                             Field(votable, name='SEMIAX2_SKY', datatype='float', unit='arcsec')]

        table.fields.extend(fieldsVO)
        table.create_arrays(len(regions.idRegions))

        for nReg in range(len(regions.idRegions)):
            row = [regions.idRegions[nReg],
                   centroids[nReg, 1]+1, centroids[nReg, 0]+1, radii[nReg], sizes[nReg]]

            if wcsMask is not None:
                row += [skyCoordsCentroids[nReg, 0], skyCoordsCentroids[nReg, 1], radiiSky[nReg]]

            row += [peaks[nReg], peaksPix[nReg, 1]+1, peaksPix[nReg, 0]+1]

            if wcsMask is not None:
                row += [skyCoordsPeak[nReg, 0], skyCoordsPeak[nReg, 1]]

            row.append(fluxes[nReg])

            if ellipse is True:
                row += [ellipseFit[nReg, 1]+1, ellipseFit[nReg, 0]+1, ellipseFit[nReg, 2], ellipseFit[nReg, 3],
                        ellipseAngle[nReg], ellipseFit[nReg, 4]]

                if wcsMask is not None:
                    row += [ellipseCentroidWCS[nReg, 0], ellipseCentroidWCS[nReg, 1], ellipseAxisWCS[nReg, 0],
                            ellipseAxisWCS[nReg, 1]]

            table.array[nReg] = tuple(row)

        votable.to_xml(voTableCat)
        logger.write('VOTable catalogue: %s' % voTableCat)

    else:
    # If the vo module is not available, uses PyFits to create a Fits table

        logger.write('Creating FITS table ... ', newLine=True)

        fields = [pf.Column(name='ID', format='J', array=regions.idRegions),
                  pf.Column(name='X_CENTROID', format='E', unit='pix', array=centroids[:, 1]+1),
                  pf.Column(name='Y_CENTROID', format='E', unit='pix', array=centroids[:, 0]+1),
                  pf.Column(name='RADIUS_PIX', format='E', unit='pix', array=radii),
                  pf.Column(name='SIZE', format='E', unit='pix', array=sizes)]

        if wcsMask is not None:
            fields += [pf.Column(name='RA_CENTROID', format='D', unit='deg',
                                 array=skyCoordsCentroids[:, 0]),
                       pf.Column(name='DEC_CENTROID', format='D', unit='deg',
                                 array=skyCoordsCentroids[:, 1]),
                       pf.Column(name='RADIUS_SKY', format='D', unit='arcsec',
                                 array=radiiSky)]

        fields += [pf.Column(name='PEAK', format='E', unit='CPS', array=peaks),
                   pf.Column(name='X_PEAK', format='E', unit='pix', array=peaksPix[:, 1]+1),
                   pf.Column(name='Y_PEAK', format='E', unit='pix', array=peaksPix[:, 0]+1)]

        if wcsMask is not None:
            fields += [pf.Column(name='RA_PEAK', format='D', unit='deg',
                                 array=skyCoordsPeak[:, 0]),
                       pf.Column(name='DEC_PEAK', format='D', unit='deg',
                                 array=skyCoordsPeak[:, 1])]

        fields += [pf.Column(name='FLUX', format='E', unit='CPS', array=fluxes)]

        if ellipse is True:
            fields += [pf.Column(name='X_ELLIPSE', format='E', unit='pix',
                                 array=ellipseFit[:, 1]+1),
                       pf.Column(name='Y_ELLIPSE', format='E', unit='pix',
                                 array=ellipseFit[:, 0]+1),
                       pf.Column(name='SEMIAX1_PIX', format='E', unit='pix',
                                 array=ellipseFit[:, 2]),
                       pf.Column(name='SEMIAX2_PIX', format='E', unit='pix',
                                 array=ellipseFit[:, 3]),
                       pf.Column(name='ELL_ANGLE', format='E', unit='deg',
                                 array=ellipseAngle),
                       pf.Column(name='ELL_ERROR', format='L', array=ellipseFit[:, 4])]

            if wcsMask is not None:
                fields += [pf.Column(name='RA_ELLIPSE', format='D', unit='deg',
                                     array=ellipseCentroidWCS[:, 0]),
                           pf.Column(name='DEC_ELLIPSE', format='D', unit='deg',
                                     array=ellipseCentroidWCS[:, 1]),
                           pf.Column(name='SEMIAX1_SKY', format='E', unit='arcsec',
                                     array=ellipseAxisWCS[:, 0]),
                           pf.Column(name='SEMIAX2_SKY', format='E', unit='arcsec',
                                     array=ellipseAxisWCS[:, 1])]

        tbhdu = pf.new_table(fields)
        hdu = pf.PrimaryHDU()
        thdulist = pf.HDUList([hdu, tbhdu])

        fitsTableCat = os.path.splitext(voTableCat)[0] + '_Cat.fits'
        if os.path.exists(fitsTableCat): os.remove(fitsTableCat)
        thdulist.writeto(fitsTableCat)
        logger.write('FITS catalogue: %s' % fitsTableCat)

    # Creates the DS9 regions file using the centroid and radius information for each region.
    logger.write('Creating DS9 regions files ... ', newLine=True)

    # For testing purposes we write a file in X, Y.
    ds9XYFile = os.path.splitext(ds9RegsFile)[0] + '_XY.reg'
    unitRegCentroid = open(ds9XYFile, 'w')
    print >>unitRegCentroid, 'global color = green'

    ii = 0
    for yy, xx in centroids:
        regID = regions.idRegions[ii]
        rad = radii[ii]
        if rad > 0.0:
            print >>unitRegCentroid, 'image; circle ' + str(xx+1) + ' ' + str(yy+1) + ' ' + \
                str(rad) + '# text = {%d}' % regID
            # print >>unitRegCentroid, 'image; x point ' + str(xx+1) + ' ' + str(yy+1) + ' # color=green'
        ii += 1

    unitRegCentroid.close()

    # Write RA, Dec instead:
    if wcsMask is not None:
        ds9WCSFile = os.path.splitext(ds9RegsFile)[0] + '_WCS.reg'
        unitRegCentroid = open(ds9WCSFile, 'w')

        print >>unitRegCentroid, 'global color = blue'
        print >>unitRegCentroid, 'fk5'

        ii = 0
        for ra, dec in skyCoordsCentroids:
            regID = regions.idRegions[ii]
            rad = radii[ii]
            if rad > 0.0:
                print >> unitRegCentroid, 'circle(' + str(ra) + ',' + str(dec) + ',' + str(rad) + '")' + \
                    ' # text = {%d}' % regID
                ii += 1

        unitRegCentroid.close()

    if ellipse:
        ds9EllFile = os.path.splitext(ds9RegsFile)[0] + '_Ell.reg'
        unitRegEllipse = open(ds9EllFile, 'w')

        print >>unitRegEllipse, 'global color = red'

        for ii in range(regions.nRegs):
            regID = regions.idRegions[ii]
            print >>unitRegEllipse, 'image; ellipse ' + str(ellipseFit[ii, 1]+1) + ' ' + \
                str(ellipseFit[ii, 0]+1) + ' ' + str(ellipseFit[ii, 2]) + ' ' + str(ellipseFit[ii, 3]) + \
                ' ' + str(ellipseAngle[ii]) + '# text = {%d}' % regID
            # print >>unitRegEllipse, 'image; x point ' + str(xx+1) + ' ' + str(yy+1) + ' # color=green'
            ii += 1

        unitRegEllipse.close()

    # Writes a file with the peak positions
    if peaksFile is not None:

        def dec2dms(dd):
            sign = -1. if dd < 0 else 1.  # Remember to deal with negative DEC
            mnt, sec = divmod(abs(dd) * 3600, 60)
            deg, mnt = divmod(mnt, 60)
            return deg * sign, mnt, sec

        logger.write('Saving peak positions ... ', newLine=True)

        peaksUnit = open(peaksFile, 'w')
        nn = 0
        for yy, xx in peaksPix:
            if wcsMask is not None:
                ra = skyCoordsPeak[nn, 0]
                dec = skyCoordsPeak[nn, 1]
                raHH, raMM, raSS = dec2dms(ra/15.)
                decDD, decMM, decSS = dec2dms(dec)
                row = '%d, %.2f %.2f, %02d %02d %05.2f, %+02d %02d %05.2f' % \
                      (nn, xx, yy, raHH, raMM, raSS, decDD, decMM, decSS)
            else:
                raiseError('No WCS information. Cannot continue.')

            print >>peaksUnit, row

            nn += 1

        peaksUnit.close()

    if plot:
        logger.write('Plotting regions ... ', newLine=True)

        try:
            from plotDS9 import plotDS9
            import ds9
        except:
            raiseWarning('No pyDS9 module find. Skipping plot.')
            return

        d = ds9.ds9()
        d.set('frame delete all')
        plotDS9(image, frame=1)
        d.set('scale mode minmax')
        d.set('regions load ' + os.path.realpath(ds9XYFile))
        plotDS9(mask, frame=2)
        d.set('scale log')
        d.set('scale mode minmax')
        d.set('regions load ' + os.path.realpath(ds9XYFile))
        if ellipse:
            plotDS9(image, frame=3)
            d.set('scale mode minmax')
            d.set('regions load ' + os.path.realpath(ds9EllFile))
            plotDS9(mask, frame=4)
            d.set('scale log')
            d.set('scale mode minmax')
            d.set('regions load ' + os.path.realpath(ds9EllFile))

    return


# Parses the command line options and call the main createCat routine
# if __name__ == '__main__':

#     from optparse import OptionParser

#     parser = OptionParser()

#     print

#     usage = 'usage: %prog [options] image mask [output] [ds9RegsFile]'
#     parser = OptionParser(usage=usage)

#     parser.add_option('-p', '--plot', dest='plot', default=False,
#         help="plots the .reg file in DS9", action="store_true")

#     parser.add_option('-e', '--ellipse', dest='ellipse', default=False,
#         help='fits ellipses to the regions', action='store_true')

#     (options, args) = parser.parse_args()

#     if len(args) == 0:
#         parser.error('incorrect number of arguments\n')
#     elif len(args) == 1:
#         image = args[0]
#         mask = os.path.splitext(image)[0] + '_Msk.fits'
#         print 'Assuming Mask file: {0}'.format(mask)
#         output = os.path.splitext(image)[0] + '.vot'
#         print 'Assuming Output file: {0}'.format(output)
#         ds9RegsFile = os.path.splitext(image)[0] + '.reg'
#         print 'Assuming DS9 regions file: {0}\n'.format(ds9RegsFile)
#     elif len(args) == 2:
#         image = args[0]
#         mask = args[1]
#         output = os.path.splitext(image)[0] + '.vot'
#         print 'Assuming Output: {0}'.format(output)
#         ds9RegsFile = os.path.splitext(image)[0] + '.reg'
#         print 'Assuming DS9 regions file: {0}\n'.format(ds9RegsFile)
#     elif len(args) == 3:
#         image = args[0]
#         mask = args[1]
#         output = args[2]
#         ds9RegsFile = os.path.splitext(image)[0] + '.reg'
#         print 'Assuming DS9 regions file: {0}\n'.format(ds9RegsFile)
#     else:
#         image = args[0]
#         mask = args[1]
#         output = args[2]
#         ds9RegsFile = args[3]

#     if False in map(os.path.exists, [image, mask]):
#         parser.error('Either the image or the mask does not exist\n')

#     createCatalog(image, mask, output, ds9RegsFile, plot=options.plot, ellipse=options.ellipse)
#     print
