#!/usr/bin/env python
# encoding: utf-8
"""
calcBackground.py

Created by José Ramón Sánchez-Gallego on 2011-02-17.
Copyright (c) 2011. All rights reserved.

"""

from plotDS9 import plotDS9
import sys, ds9
import pyfits as pf
import numpy as np

def calcBack(data, sigma=3.):
    # print
    # print "  #       Mean          Std           Dif."
    tmpData = np.array(data)
    mean = np.mean(tmpData); std = np.std(tmpData)
    step = 0
    while step < 1000:
        oldMean = mean; oldStd = std
        minValue = oldMean - sigma * oldStd
        maxValue = oldMean + sigma * oldStd
        tmpData = tmpData[np.where(tmpData >= minValue)]
        tmpData = tmpData[np.where(tmpData <= maxValue)]
        mean = np.mean(tmpData); std = np.std(tmpData)
        dif = np.abs(((mean-oldMean)/oldMean)*100.0)
        # print '%3d   %11.4E   %11.4E   %11.4E' % (step+1, mean, std, dif)
        if dif < 0.02: return mean
        step += 1 
    return mean


def calcBackground(image, fancyPrint):
    
    import tempfile

    fancyPrint('', newLine=True)
    
    sys.stdout.write('Displaying image and background region ... ')
    sys.stdout.flush()
    
    d = ds9.ds9()
    d.set('frame delete all')
    plotDS9(image, frame=1)
    
    xSize,ySize = map(int,d.get('fits size').split())

    fileTmpReg = tempfile.NamedTemporaryFile()
    print >>fileTmpReg, 'image; box %d %d %d %d # color = green' % (xSize/2, ySize/2, xSize/30, ySize/30)
    fileTmpReg.flush()

    d.set('regions system image')
    d.set('regions format ds9')
    d.set('regions load ' + fileTmpReg.name)
    
    fileTmpReg.close()

    print 'Done'
    
    raw_input('Select the region which will be used to estimate the background ... ')
    
    d.set('regions select all')
    regionOut = d.get('regions')
    pos0 = regionOut.find('box(') + 4
    xCenter, yCenter, xSize, ySize = map(float, regionOut[pos0:].split(',')[0:-1])
    x0 = xCenter - xSize/2.; xf = xCenter + xSize/2.
    y0 = yCenter - ySize/2.; yf = yCenter + ySize/2.
    
    # sys.stdout.write('Calculating background level ... ')
    # sys.stdout.flush()
    data = pf.getdata(image)
    print data.shape
    dataBack = data[y0-1:yf-1,x0-1:xf-1]
    backLevel = calcBack(dataBack)
    # print '%11.5e' % backLevel
    
    fancyPrint('')
    
    return backLevel
