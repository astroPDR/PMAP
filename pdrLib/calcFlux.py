#!/usr/bin/env python
# encoding: utf-8
"""

calcFlux.py

Created by José Ramón Sánchez-Gallego on 2011-08-03.
Copyright (c) 2011. All rights reserved.

"""


import sys, os, re
import tempfile
from pyraf import iraf
import numpy as np
import pyfits #to check for EXPTIME keyword

"""
Given a image name, this script calculates the flux within a circular aperture.
The script uses pyraf and does not performs any kind of recentering or background
estimation.
"""
def calcFlux(imageName, xx, yy, rr):
    
    coordsFN = tempfile.NamedTemporaryFile(delete=False)
    print >>coordsFN, xx, yy
    coordsFN.close()
    
    digiphot = iraf.digiphot(_doprint=0)
    apphot = iraf.apphot(_doprint=0)
    qphot = iraf.qphot
    pdump = iraf.pdump
    
    qphot.coords = coordsFN.name
    
    #if (iraf.hedit(image=imageName,fields='EXPTIME',value='.')):
    #^this line prints the value of EXPTIME which is not needed
    #need to check for EXPTIME or qphot will raise an error
    try:
      pyfits.open(imageName)[0].header['EXPTIME']
      qphot.exposure = 'EXPTIME'
    except:
      pass
    
    qphot.airmass = ''
    
    qphot.interactive = False
    qphot.radplots = False
    qphot.verbose = False
    
    annulus = rr
    dannulus = 1.0
    
    qphot.output = 'qPhotOutput.tmp'
    if os.path.exists(qphot.output): os.remove(qphot.output)
    qphot(imageName, 1.0, annulus, dannulus, annulus, Stdout=1)
    
    areas = np.array(map(float, pdump(qphot.output, 'AREA', 'yes', Stdout=1)[0].split()))
    fluxes = np.array(map(float, pdump(qphot.output, 'SUM', 'yes', Stdout=1)[0].split()))
    
    if os.path.exists(qphot.output): os.remove(qphot.output)
    
    return [np.float64(areas[0]), np.float64(fluxes[0])]



"""
This routine calculates the aperture photometry for the series of positions given in coords.
The photometry is done for each aperture from rMin to rMax each step pixels. Recentering
is allowed and controlled through centerBox.
"""
def calcFluxMulti(imageName, coords, rMin, rMax, step=1, centerBox=5.0):
    
    nCoords = coords.shape[0]
    
    results = []
    
    coordsFN = tempfile.NamedTemporaryFile(delete=False)
    for xx, yy in coords: print >>coordsFN, xx, yy
    coordsFN.close()
    
    digiphot = iraf.digiphot(_doprint=0)
    apphot = iraf.apphot(_doprint=0)
    qphot = iraf.qphot
    pdump = iraf.pdump
    
    qphot.coords = coordsFN.name
    
    #if (iraf.hedit(image=imageName,fields='EXPTIME',value='.')):
    #^this line prints the value of EXPTIME which is not needed
    #need to check for EXPTIME or qphot will raise an error
    try:
      pyfits.open(imageName)[0].header['EXPTIME']
      qphot.exposure = 'EXPTIME'
    except:
      pass
    
    qphot.airmass = ''
    
    qphot.interactive = False
    qphot.radplots = False
    qphot.verbose = False
    
    apertures = ','.join([str(nn) for nn in range(int(rMin), int(rMax+1))])
    annulus = rMax
    dannulus = rMax + 1
    
    qphot.output = 'qPhotOutput.tmp'
    if os.path.exists(qphot.output): os.remove(qphot.output)
    qphot(imageName, centerBox, annulus, dannulus, apertures, Stdout=1)
    
    for nn in range(nCoords):
    
        areas = np.array(map(float, pdump(qphot.output, 'AREA', 'yes', Stdout=1)[nn].split()))
        fluxes = np.array(map(float, pdump(qphot.output, 'SUM', 'yes', Stdout=1)[nn].split()))
        newXX = np.array(map(float, pdump(qphot.output, 'XCENTER', 'yes', Stdout=1)[nn].split()))[0]
        newYY = np.array(map(float, pdump(qphot.output, 'YCENTER', 'yes', Stdout=1)[nn].split()))[0]
        results.append([nn, coords[nn][0], coords[nn][1], newXX, newYY, areas, fluxes])
        
    if os.path.exists(qphot.output): os.remove(qphot.output)
    
    return results
    
    

    
