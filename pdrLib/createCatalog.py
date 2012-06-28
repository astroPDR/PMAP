#!/usr/bin/env python
# encoding: utf-8
"""
createCatalog.py

Created by José Ramón Sánchez-Gallego on 2011-03-25.
Copyright (c) 2011. All rights reserved.

This script creates a VOTable catalog with the properties of each region within
the selected mask.

"""

#2011-10-17: updated to deal with negative DEC correctly

import numpy as np
import sys, os
import pyfits as pf
import pywcs as pw
from fitellipse import fitellipse
from Regions import *
import atpy

def flush(text):
    sys.stdout.write(text)
    sys.stdout.flush()


# Main routine which read the image and mask files, calculates various parameters and
# creates the VOTable catalogue and the DS9 regions files.
def createCatalog(image, mask, output, ds9RegsFile, peaksFile=None, ellipse=False, ellMode='linear', plot=False):
    
# Reads the header, data and loads the astrometric solution
    hduImage = pf.open(image)
    hduMask = pf.open(mask)
    try:
        wcsMask = pw.WCS(hduMask[0].header)
    except:
        wcsMask == None
        print 'No WCS information. The VOTable file will not include those fields.\n'
    
    flush('Gathering mask data ... ')
    regions = RegionSet(hduMask[0].data, image=hduImage[0].data)
    print 'done'
    
    flush('Preparing data ... ')
    centroids = np.array([regions.Regions[idReg].getCentroid()[0] 
        for idReg in regions.idRegions])
    radii = [regions.Regions[idReg].getCentroid()[1] for idReg in regions.idRegions]
    sizes = [regions.Regions[idReg].where[0].size for idReg in regions.idRegions]
    peaks = [regions.Regions[idReg].peak for idReg in regions.idRegions]
    peaksPix = np.array([regions.Regions[idReg].peakWhere 
        for idReg in regions.idRegions])+1
    fluxes = [regions.Regions[idReg].getFlux() for idReg in regions.idRegions]
    
    if wcsMask != None:
        skyCoordsCentroids = wcsMask.wcs_pix2sky(np.fliplr(centroids), 0)
        skyCoordsPeak = wcsMask.wcs_pix2sky(np.fliplr(peaksPix), 0)
        
        try:
            scale = np.abs(hduImage[0].header['CD1_1']) * 3600.
        except:
            scale = np.abs(hduImage[0].header['CDELT1']) * 3600.
        radiiSky = radii * scale
    
    print 'done'
    
# Calculates an ellipse fit for the pixels in the region
    if ellipse:
        flush('Fitting ellipse parameters ... ')
        def getEllParams(data, id):
            try:
                return list(fitellipse(data, ellMode)) + [0]
            except RuntimeError:
                return [regions.Regions[id].getCentroid()[0][0], regions.Regions[id].getCentroid()[0][1],
                    regions.Regions[idReg].getCentroid()[1], regions.Regions[idReg].getCentroid()[1], 0.0, 1]
        ellipseFit = np.array([getEllParams(regions.Regions[idReg].pixels, idReg) for idReg in regions.idRegions])        
        ellipseCentroidWCS = np.array(wcsMask.wcs_pix2sky(ellipseFit[:,[1,0]], 0))
        ellipseAxisWCS = np.array(ellipseFit[:,[2,3]] * scale)
        ellipseAngle = np.array(-ellipseFit[:,4] * 180. / np.pi + 90.)
        print 'done'
    
    flush('Creating VOTable ... ')
    
    if os.path.exists(output): os.remove(output)
    
    t = atpy.Table()
    t.add_column('ID', regions.idRegions, dtype=np.int16)
    t.add_column('X_CENTROID', centroids[:,1] + 1., unit='pixels', dtype=np.float32)
    t.add_column('Y_CENTROID', centroids[:,0] + 1., unit='pixels', dtype=np.float32)
    t.add_column('RADIUS_PIX', radii, unit='pixels', dtype=np.float32)
    t.add_column('SIZE', sizes, unit='pixels', dtype=np.float32)
    if wcsMask != None:
        t.add_column('RA_CENTROID', skyCoordsCentroids[:,0], unit='degrees', dtype=np.float32)
        t.add_column('DEC_CENTROID', skyCoordsCentroids[:,1], unit='degrees', dtype=np.float32)
        t.add_column('RADIUS_SKY', radiiSky, unit='arcsec', dtype=np.float32)
    t.add_column('PEAK', peaks, unit='counts', dtype=np.float32)
    t.add_column('X_PEAK', peaksPix[:,1] + 1., unit='pixels', dtype=np.float32)
    t.add_column('Y_PEAK', peaksPix[:,0] + 1., unit='pixels', dtype=np.float32)
    if wcsMask != None:
        t.add_column('RA_PEAK', skyCoordsPeak[:,0], unit='degrees', dtype=np.float32)
        t.add_column('DEC_PEAK', skyCoordsPeak[:,1], unit='degrees', dtype=np.float32)
    t.add_column('FLUX', fluxes, unit='counts', dtype=np.float32)
    if ellipse:
        t.add_column('X_ELLIPSE', ellipseFit[:,1] + 1., unit='pixels', dtype=np.float32)
        t.add_column('Y_ELLIPSE', ellipseFit[:,0] + 1., unit='pixels', dtype=np.float32)
        t.add_column('SEMIAX1_PIX', ellipseFit[:,2], unit='pixels', dtype=np.float32)
        t.add_column('SEMIAX2_PIX', ellipseFit[:,3], unit='pixels', dtype=np.float32)
        t.add_column('RA_ELLIPSE', ellipseCentroidWCS[:,0], unit='degrees', dtype=np.float32)
        t.add_column('DEC_ELLIPSE', ellipseCentroidWCS[:,1], unit='degrees', dtype=np.float32)
        t.add_column('SEMIAX1_SKY', ellipseAxisWCS[:,0], unit='arcsec', dtype=np.float32)
        t.add_column('SEMIAX2_SKY', ellipseAxisWCS[:,1], unit='arcsec', dtype=np.float32)
        t.add_column('ELL_ANGLE', ellipseAngle, unit='degrees', description='Given as for the DS9 regions', dtype=np.int)    
        t.add_column('ELL_ERROR', ellipseFit[:,4], unit='Bool', description='0: No error, 1: Error fitting ellipse', dtype=np.int)    
        
    t.write(output, verbose=False) # Saves it
    
    print 'done'
    
# Creates the DS9 regions file using the centroid and radius information for each region.
    flush('Creating DS9 regions files ... ')
    
    # For testing purposes we write a file in x, y and ra, dec
    unitRegCentroid = open('xy_' + ds9RegsFile, 'w')
    print >>unitRegCentroid, 'global color = green'
    ii = 0
    for yy,xx in centroids:
        rad = radii[ii]
        if rad > 0.0:
            print >>unitRegCentroid, 'image; circle ' + str(xx+1) + ' ' + str(yy+1) + ' ' + str(rad)
            # print >>unitRegCentroid, 'image; x point ' + str(xx+1) + ' ' + str(yy+1) + ' # color=green'
        ii += 1
    unitRegCentroid.close()
    # Write ra, dec instead:
    unitRegCentroid = open(ds9RegsFile, 'w')
    print >>unitRegCentroid, 'global color = blue'
    print >>unitRegCentroid, 'fk5'
    ii = 0
    for ra,dec in skyCoordsCentroids:
          rad = radii[ii]
          if rad > 0.0:
                print >> unitRegCentroid, 'circle(' + str(ra) + ',' + str(dec) + ',' + str(rad) + '")'
                ii += 1
    unitRegCentroid.close()
    
    if ellipse:
        ds9EllFile = os.path.splitext(ds9RegsFile)[0] + '_Ell.reg'
        unitRegEllipse = open(ds9EllFile, 'w')
        print >>unitRegEllipse, 'global color = green'
        for ii in range(regions.nRegs):
            print >>unitRegEllipse, 'image; ellipse ' + str(ellipseFit[ii,1]+1) + ' ' + \
                str(ellipseFit[ii,0]+1) + ' ' + str(ellipseFit[ii,2]) + ' ' + str(ellipseFit[ii,3]) + \
                ' ' + str(ellipseAngle[ii])
            # print >>unitRegEllipse, 'image; x point ' + str(xx+1) + ' ' + str(yy+1) + ' # color=green'
            ii += 1
        unitRegEllipse.close()
        
    print 'done'
    
    if peaksFile != None:
        
        def dec2dms(dd):
            sign = -1. if dd < 0 else 1. #Remember to deal with negative DEC
            mnt,sec = divmod(abs(dd)*3600,60)
            deg,mnt = divmod(mnt,60)
            return deg*sign,mnt,sec
        
        flush('Saving peak positions ... ')
        
        peaksUnit = open(peaksFile, 'w')
        nn = 1
        for ra,dec in skyCoordsPeak:
            raHH, raMM, raSS = dec2dms(ra/15.)
            decDD, decMM, decSS = dec2dms(dec)
            row = '%d, %02d %02d %05.2f, %+02d %02d %05.2f' % (nn, raHH, raMM, raSS, decDD, decMM, decSS)
            print >>peaksUnit, row
            nn += 1
        peaksUnit.close()
        
        print 'done'
        
    
    if plot:
        flush('Plotting regions ... ')
        from plotDS9 import plotDS9
        import ds9
        d = ds9.ds9()
        plotDS9(image, frame=1)
        d.set('regions load ' + os.path.realpath(ds9RegsFile))
        plotDS9(mask, frame=2)
        d.set('scale log')
        d.set('scale mode minmax')
        d.set('regions load ' + os.path.realpath(ds9RegsFile))
        if ellipse:
            plotDS9(image, frame=3)
            d.set('regions load ' + os.path.realpath(ds9EllFile))
            plotDS9(mask, frame=4)
            d.set('scale log')
            d.set('scale mode minmax')
            d.set('regions load ' + os.path.realpath(ds9EllFile))
        print 'done'
        
    return
    
  
# Parses the command line options and call the main createCat routine  
if __name__ == '__main__':

    from optparse import OptionParser

    parser = OptionParser()
    
    print
    
    usage = 'usage: %prog [options] image mask [output] [ds9RegsFile]'
    parser = OptionParser(usage=usage)
    
    parser.add_option('-p', '--plot', dest='plot', default=False,
        help="plots the .reg file in DS9", action="store_true")

    parser.add_option('-e', '--ellipse', dest='ellipse', default=False,
        help='fits ellipses to the regions', action='store_true')
    
    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error('incorrect number of arguments\n')
    elif len(args) == 1:
        image = args[0]
        mask = os.path.splitext(image)[0] + '_Msk.fits'
        print 'Assuming Mask file: {0}'.format(mask)
        output = os.path.splitext(image)[0] + '.vot'
        print 'Assuming Output file: {0}'.format(output)
        ds9RegsFile = os.path.splitext(image)[0] + '.reg'
        print 'Assuming DS9 regions file: {0}\n'.format(ds9RegsFile)
    elif len(args) == 2:
        image = args[0]
        mask = args[1]
        output = os.path.splitext(image)[0] + '.vot'
        print 'Assuming Output: {0}'.format(output)
        ds9RegsFile = os.path.splitext(image)[0] + '.reg'
        print 'Assuming DS9 regions file: {0}\n'.format(ds9RegsFile)
    elif len(args) == 3:
        image = args[0]
        mask = args[1]
        output = args[2]
        ds9RegsFile = os.path.splitext(image)[0] + '.reg'
        print 'Assuming DS9 regions file: {0}\n'.format(ds9RegsFile)
    else:
        image = args[0]
        mask = args[1]
        output = args[2]
        ds9RegsFile = args[3]
    
    if False in map(os.path.exists, [image, mask]):
        parser.error('Either the image or the mask does not exist\n')
    
    createCatalog(image, mask, output, ds9RegsFile, plot=options.plot, ellipse=options.ellipse) 
    print
