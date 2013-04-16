#!/usr/bin/env python
# encoding: utf-8
"""
calcBackground.py

Created by José Sánchez-Gallego on 2013-04-15.
Copyright (c) 2013. All rights reserved.

This routine calculates the background level of an image using the
following process:
a) A temporary background level of the image is calculated by removing
   all pixels above 3-sigma and then calculating the mean value.
b) CLFind is run on the image using this background noise and generic
   values for the contour levels.
c) Using the mask generated with CLFind, the residual image is calculated.
d) Steps a) is repeated and the background level determined that way is
   returned.

"""

import sys
from clfind2d import clfind2d
import pyfits as pf
import numpy as np
import os


def calcBack(data, sigma=3.):

    # print
    # print "  #       Mean          Std           Dif."

    tmpData = data.flatten()[np.isfinite(data.flatten())]

    mean = np.mean(tmpData)
    std = np.std(tmpData)
    step = 0
    while step < 1000:
        oldMean = mean
        oldStd = std
        minValue = oldMean - sigma * oldStd
        maxValue = oldMean + sigma * oldStd
        tmpData = tmpData[np.where(tmpData >= minValue)]
        tmpData = tmpData[np.where(tmpData <= maxValue)]
        mean = np.mean(tmpData)
        std = np.std(tmpData)
        dif = np.abs(((mean-oldMean)/oldMean)*100.0)
        if dif < 0.02: return mean
        step += 1
    return mean


def calcBackground(image, logger, verbose=True):

    import tempfile

    print
    print 'Automatic estimation of the background level'

    imageData = pf.open(image)

    tmpBackgroundLevel = calcBack(imageData[0].data)
    print 'Step 1) Background level = %.4E' % tmpBackgroundLevel

    sys.stdout.write('Step 2) Running CLFind to obtain mask ... ')
    sys.stdout.flush()

    if verbose is True:
        print

    tmpMask = tempfile.NamedTemporaryFile(dir='./', delete=False)
    tmpLevels = np.linspace(10.0 * tmpBackgroundLevel, np.max(imageData[0].data), 20)

    clfind2d(image, tmpMask.name, tmpLevels, log=False, nPixMin=20, verbose=verbose)

    if verbose is False:
        print 'done'

    print 'Step 3) Calculating residuals'
    dataMask = pf.getdata(tmpMask.name)
    mask = np.where(dataMask > 0)
    imageData[0].data[mask] = np.nan

    backLevel = calcBack(imageData[0].data)
    print 'Step 4) Background level of the residuals = %.4E' % backLevel

    if os.path.exists(tmpMask.name):
        os.remove(tmpMask.name)

    del imageData, dataMask

    return backLevel
