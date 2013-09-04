#!/usr/bin/env python
# encoding: utf-8
"""
fixTHINGS.py

Created by José Ramón Sánchez-Gallego on 2010-06-14.
Copyright (c) 2010. All rights reserved.

Fix the header of THINGS files to remove the references to axis 3 and 4.
It solves problems when using those images with pyfits and pywcs.

"""

import sys
import os
from astropy.io import fits as pf


def fixTHINGS(imageIn, imageOut):

    print
    sys.stdout.write('Fixing file %s ... ' % imageIn)
    sys.stdout.flush()

    if imageOut != imageIn:
        hdu = pf.open(imageIn)
    else:
        hdu = pf.open(imageIn, mode='update')

    dataNew = hdu[0].data[0, 0, :, :]

    del hdu[0].header['CTYPE3']
    del hdu[0].header['CDELT3']
    del hdu[0].header['CRVAL3']
    del hdu[0].header['CRPIX3']
    del hdu[0].header['CROTA3']

    del hdu[0].header['CTYPE4']
    del hdu[0].header['CDELT4']
    del hdu[0].header['CRVAL4']
    del hdu[0].header['CRPIX4']
    del hdu[0].header['CROTA4']

    if imageOut != imageIn:
        if os.path.exists(imageOut):
            os.remove(imageOut)
        pf.writeto(imageOut, dataNew, hdu[0].header)
    else:
        hdu[0].data = dataNew
        hdu.flush()

    print 'Done'

    print
    return


if __name__ == '__main__':

    imageIn = sys.argv[1]
    imageOut = sys.argv[2]

    fixTHINGS(imageIn, imageOut)
