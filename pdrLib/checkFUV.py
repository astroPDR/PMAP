#!/usr/bin/env python
# encoding: utf-8
"""
checkFUV.py

Created by José Ramón Sánchez-Gallego on 2011-03-08.
Copyright (c) 2011. All rights reserved.

This script localized regions that are "touching" and checks them to make sure
that they are two different regions. To do that:
(1) Finds the peak position for each test region
(2) Localizes the pixels joining the peak positions with a straight line
(3) Calculates the intensity profile across that line
(4) If certain conditions are not fulfilled, both test regions are joined
in a single region

"""

import numpy as np
import sys, os
import pyfits as pf
import pywcs as pw
import Regions

def flush(text):
    sys.stdout.write(text)
    sys.stdout.flush()
    return


def calcAngDist(coords1, coords2):
# Taken from astCoords
    
    RADeg1 = coords1[0]
    RADeg2 = coords2[0]
    decDeg1 = coords1[1]
    decDeg2 = coords2[1]
    
    cRA=np.radians(RADeg1) 
    cDec=np.radians(decDeg1) 
    
    gRA=np.radians(RADeg2) 
    gDec=np.radians(decDeg2) 
    
    dRA=cRA-gRA 
    dDec=gDec-cDec 
    cosC=(np.sin(gDec)*np.sin(cDec))+(np.cos(gDec)*np.cos(cDec)*np.cos(gRA-cRA)) 
    x=(np.cos(cDec)*np.sin(gRA-cRA))/cosC 
    y=((np.cos(gDec)*np.sin(cDec))-(np.sin(gDec)*np.cos(cDec)*np.cos(gRA-cRA)))/cosC 
    r=np.degrees(np.sqrt(x*x+y*y))
    
    return r


def checkFUV(options, fancyPrint):
    
    fuvImage = options['fuvImage']
    fuvMask = options['fuvMaskFile']
    
    fancyPrint('Checking FUV regions ... ', newLine=True, log=True)
    
    fancyPrint('Gathering regions ... ', flush=True)
    
    hduMask = pf.open(fuvMask)
    hduImage = pf.open(fuvImage)
    nRegs = np.max(hduMask[0].data)
    regions = Regions.RegionSet(hduMask[0].data, image=hduImage[0].data)
    fancyPrint('done')
    
    # fancyPrint('Getting adjacent regions ... ', flush=True)
    adjacentRegions = regions.getAdjacentRegs()
    # fancyPrint('done')
    
    fancyPrint('Processing adjacent regions ... ', flush=True)
    
    while adjacentRegions != []:
        
        pair = adjacentRegions[0]
        fluxA = regions.Regions[pair[0]].getFlux()
        fluxB = regions.Regions[pair[1]].getFlux()
        peakA = regions.Regions[pair[0]].peak
        peakB = regions.Regions[pair[1]].peak
        
        if options['fluxContrast'] != False:
            contrast = float(np.max([fluxA, fluxB]) / np.min([fluxA, fluxB]))
            minRegion = pair[np.argmin([fluxA, fluxB])]
            maxRegion = pair[np.argmax([fluxA, fluxB])]
            if options['fluxContrast'] == True or contrast > options['fluxContrast']:
                if options['joinRegions']:
                    regions.joinRegions(maxRegion, minRegion)
                    fancyPrint('Adjacent regions {0} and {1} joined. Flux contrast={2:.1f}'.format(maxRegion, minRegion, contrast), noPrint=True, log=True)
                else:
                    regions.deleteReg(minRegion)
                    fancyPrint('Adjacent regions {0} deleted. Flux contrast={1:.1f}'.format(minRegion, contrast), noPrint=True, log=True)
                adjacentRegions = regions.getAdjacentRegs()
                continue            
                
        if options['peakContrast'] != False:
            contrast = float(np.max([peakA, peakB]) / np.min([peakA, peakB]))
            minRegion = pair[np.argmin([peakA, peakB])]
            maxRegion = pair[np.argmax([peakA, peakB])]
            if options['peakContrast'] == True or contrast > options['peakContrast']:
                if options['joinRegions']:
                    regions.joinRegions(maxRegion, minRegion)
                    fancyPrint('Adjacent regions {0} and {1} joined. Peak contrast={2:.1f}'.format(maxRegion,  minRegion, contrast), noPrint=True, log=True)
                else:
                    regions.deleteReg(minRegion)
                    fancyPrint('Adjacent regions {0} deleted. Peak contrast={1:.1f}'.format(minRegion, contrast), noPrint=True, log=True)
                adjacentRegions = regions.getAdjacentRegs()
                continue
        del adjacentRegions[0]
    fancyPrint('done')


    # fancyPrint('Getting close regions ... ', flush=True)
    if options['scale'] != None:
        scale = options['scale']
    else:
        try:
            try:
                scaleDeg = np.abs(hduImage[0].header['CD1_1'])
            except:
                scaleDeg = np.abs(hduImage[0].header['CDELT1'])
        except:
            fancyPrint('failed')
            fancyPrint('Pixel scale cannot be calculated. Use the scale parameter.', newLine=True, log=True)
            print
            sys.exit()
        try:
            scalePc = 2.0 * options['distance'] * 1e6 * np.tan(0.5 * scaleDeg * np.pi / 180.)
        except:
            fancyPrint('failed')
            fancyPrint('Distance to the galaxy not defined.', newLine=True, log=True)
            print
            sys.exit()
            
    minNumPixels = options['minDistance'] / scalePc
    
    closeRegions = regions.getCloseRegs(minNumPixels, adjacent=False)
    # fancyPrint('done')
    
    fancyPrint('Processing close regions ... ', flush=True)
    
    while closeRegions != []:
        
        pair = closeRegions[0]
        
        fluxA = regions.Regions[pair[0]].getFlux()
        fluxB = regions.Regions[pair[1]].getFlux()
        peakA = regions.Regions[pair[0]].peak
        peakB = regions.Regions[pair[1]].peak
        
        if options['alwaysRemoveClose']:
            minRegion = pair[np.argmin([fluxA, fluxB])]
            regions.deleteReg(minRegion)
            closeRegions = regions.getCloseRegs(minNumPixels, adjacent=False)
            fancyPrint('Close region {0} removed'.format(minRegion), log=True, noPrint=True)
            continue
        
        if options['fluxContrast'] != False:
            contrast = float(np.max([fluxA, fluxB]) / np.min([fluxA, fluxB]))
            minRegion = pair[np.argmin([fluxA, fluxB])]
            maxRegion = pair[np.argmax([fluxA, fluxB])]
            if (options['fluxContrast'] == True) or (contrast > options['fluxContrast']):
                regions.deleteReg(minRegion)
                closeRegions = regions.getCloseRegs(minNumPixels, adjacent=False)
                fancyPrint('Close region {0} removed. Flux contrast={1:.1f}'.format(minRegion, contrast), log=True, noPrint=True)
                continue            
                
        if options['peakContrast'] != False:
            contrast = float(np.max([peakA, peakB]) / np.min([peakA, peakB]))
            minRegion = pair[np.argmin([peakA, peakB])]
            maxRegion = pair[np.argmax([peakA, peakB])]
            if (options['peakContrast'] == True) or (contrast > options['peakContrast']):
                regions.deleteReg(minRegion)
                closeRegions = regions.getCloseRegs(minNumPixels, adjacent=False)
                fancyPrint('Close region {0} removed. Peak contrast={1:.1f}'.format(minRegion, contrast), log=True, noPrint=True)
                continue
        
        del closeRegions[0]
    
    fancyPrint('done')
    
    fancyPrint('Saving new mask ... ', flush=True)
    
    hdu = pf.PrimaryHDU(regions.DataMask)
    hdu.header = hduMask[0].header
    hduList = pf.HDUList([hdu])
    if os.path.exists(options['fuvMaskFileRej']): os.remove(options['fuvMaskFileRej'])
    hduList.writeto(options['fuvMaskFileRej'], output_verify='ignore')
    
    fancyPrint('done')
    
    return


# Parses the command line options and call the main routine  
# if __name__ == '__main__':
    # checkFUV() 
