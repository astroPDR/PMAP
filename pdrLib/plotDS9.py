#!/usr/bin/env python
# encoding: utf-8
"""
plotDS9.py

Created by José Ramón Sánchez-Gallego on 2011-02-17.
Copyright (c) 2011. All rights reserved.

"""

from Error import raiseError

try:
    import ds9
except:
    raiseError('No ds9 module found.')

import os


def plotDS9(image, frame=1, zscale=True, fit=True, ext=None, adjust=True,
            minLimZero=False, keepScale=False, scale=False, delete=False):
    ddd = ds9.ds9()

    if delete: ddd.set('frame delete all')
    ddd.set('frame ' + str(frame))

    if keepScale is True:
        zMin, zMax = map(float, ddd.get('scale limits').split())

    if adjust is True:
        ddd.set('cmap heat')
        ddd.set('cmap value 1. 0.5')
    ddd.set('file ' + os.path.realpath(image))

    if fit is True: ddd.set('zoom to fit')

    if keepScale is True:
        ddd.set('scale mode zscale')
        ddd.set('scale limits %f %f' % (zMin, zMax))
        zscale = False
        scale = False

    if zscale is True:
        ddd.set('scale mode zscale')
        ddd.set('zscale contrast .25')
        ddd.set('zscale sample 600')
        ddd.set('zscale line 120')
        ddd.set('scale linear')
        if minLimZero is True:
            zMin, zMax = map(float, ddd.get('scale limits').split())
            ddd.set('scale limits %f %f' % (0.0, zMax))

    if scale is not False:
        try:
            zMin = float(scale[0])
            zMax = float(scale[1])
        except:
            print 'Problem found parsing scale values'
            return
        ddd.set('scale limits %f %f' % (zMin, zMax))

    return
