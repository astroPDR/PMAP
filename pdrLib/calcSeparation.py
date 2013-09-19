#!/usr/bin/env python
# encoding: utf-8
"""
calcSeparation.py

Created by Jose Sanchez-Gallego on 30 Aug 2013.
Based on previous code by Jonathan Heiner.
Copyright (c) 2013. All rights reserved.

Calculates the separation between two points in a galaxy taking
into account its inclination.
"""


from astropy import coordinates as coords
from astropy import units as uu
from Error import raiseError
import numpy as np
import math as m


def e_angle_rad(x1, y1, x2, y2, pa):
    # Returns angles between two vectors, taking pa into account
    # Since the output is used in a (symmetrical) ellipse equation,
    # we don't care about the quadrants
    # Assume that the first set is the center, but does it matter?
    # need to verify, but currently first set is always center

    # 1. 'center' unit vector, rotated by pa
    xc = m.cos(m.radians(pa + 90))
    yc = m.sin(m.radians(pa + 90))
    # 2. construct vector from first coord to second coord
    xp = x2 - x1
    yp = y2 - y1
    # 3. calculate angle between 1. and 2.
    lc = m.sqrt(xc ** 2. + yc ** 2.)
    lp = m.sqrt(xp ** 2. + yp ** 2.)
    # print x1, y1, x2, y2, xc, yc, xp, yp, "input values\n", lc, lp, "lc lp"
    x = (xc * xp + yc * yp) / (lc * lp)
    y = m.sqrt(1 - x ** 2.)
    # print x, y, m.atan2(y,x), m.degrees(m.atan2(y,x)), "angle between
    # vectors"
    return(m.atan2(y, x))


def createCoord(ra, dec):

    if np.isscalar(ra) is True and np.isscalar(dec) is True:
        if np.isreal(ra) is True and np.isreal(dec):
            cc = coords.ICRSCoordinates(ra, dec, unit=(uu.degree, uu.degree))
        else:
            try:
                cc = coords.ICRSCoordinates('%s %s' % (ra, dec))
            except:
                try:
                    raH, raM, raS = ra.split()
                    decD, decM, decS = dec.split()
                    cc = coords.ICRSCoordinates('%sh%sm%ss %sd%sm%ss' % (raH, raM, raS,
                                                                         decD, decM, decS))
                except:
                    raiseError('Coordinate cannot be converted to astropy.coordinates.')
    return cc


def calcSeparation(ra1, dec1, ra2, dec2, pa, incl, D_gal, wcs):

    # Creates the astropy coordinates
    c1 = createCoord(ra1, dec1)
    c2 = createCoord(ra2, dec2)

    # 1. calculate separation angle between the two coordinates, giving us r_raw
    angle = c1.separation(c2)  # In radians
    r_raw = D_gal * angle.radians

    # 2. correct for pa, i by deprojecting radius
    # For this we need x and y coordinates through the wcs header

    x1, y1 = wcs.wcs_world2pix(np.array([[c1.ra.degrees, c1.dec.degrees]]), 0)[0]
    x2, y2 = wcs.wcs_world2pix(np.array([[c2.ra.degrees, c2.dec.degrees]]), 0)[0]
    e_a = e_angle_rad(x1, y1, x2, y2, pa)  # get the angle
    esqd = 1 - (m.cos(m.radians(incl))) ** 2.  # deproject
    sqrtf = m.sqrt((1 - esqd) / (1 - esqd * m.cos(e_a) ** 2.))
    # print esqd, m.cos(e_a), sqrtf, "esqd cos(angle) sqrtf"
    # calculate separation in parsec from angle and ellipse deprojection

    return r_raw / sqrtf  # in units of D_gal (e.g. kpc)
