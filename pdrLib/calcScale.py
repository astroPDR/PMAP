#!/usr/bin/env python
# encoding: utf-8
"""
calcScale.py

Created by José Sánchez-Gallego on 26 Sep 2013.
Copyright (c) 2013. All rights reserved.

Description goes here.

Required input:
    - fileName: the filename of the image of which we want to
        get the pixel scale.

Returns:
    The average of the pixel scale in both spatial dimensions

v0.1: 26 Sep 2013 J. Sánchez-Gallego
    Initial version

"""

from . import raiseError
import numpy as np
from pdrLib import open_image


def calcScale(fileName):

    wcs, data, hdu = open_image(fileName)

    try:
        pixSize = np.array([np.abs(wcs.wcs.cd[0, 0]), np.abs(wcs.wcs.cd[1, 1])]) * 3600
    except:
        try:
            pixSize = np.abs(wcs.wcs.cdelt[0:2]) * 3600
        except:
            raiseError('Pixel size for image %s cannot be determined. ' +
                       'Please, specify it in the configuration file' % fileName)

    return pixSize.mean()
