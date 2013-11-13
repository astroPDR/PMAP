#!/usr/bin/env python
# encoding: utf-8

# Collection of functions for various parts of the PDR method 'pipeline' script pipeline.py
# Depends on astLib (0.4.0), Pyfits, NumPy & SciPy, (kapteyn), matplotlib, sextutils
# JSH 2010-02-01
#    2011-07-14
#    2012-02-14~latest edit (needed to fix a bug on Fedora with formatting)
#    2012-02-15 edits to avoid divisions by zero (should disqualify the measurement)
#                 also, fix errors when going outside the dust map
# & José Sanchez-Gallego
#    2012-06-04 math range errors for high values of 12+log(O/H)

# Fedora bug: https://bugzilla.redhat.com/show_bug.cgi?format=multiple&id=654269
#   Workaround is to 'dumb down' the formatting statements and do a type cast in advance
#   It seems to happen with 'x' and 'd' formatters
# Another workaround is to replace the 'd' with 'g', which has the desired
# effect

# pdrLib is now organized to contain José's routines as well, seperate
# from this file.

# Wishlist: classes like coordinate lists with built in options to query hms dms, or decimal etc.
#		because now we have to do all kinds of conversions first
#		can get rid of passing the header explicitely all the time

# from __future__ import with_statement #Needed for logdump's 'with' use
# in Python 2.5 only?
import pdb  # debugger
import numpy as np
from . import pf, pw
from string import strip
# import convCoords as astCoords
# for flux specifically and plot routines in general
from matplotlib import pylab as plt
import math as m
from Error import raiseError

""" Default file operations """


def read_dust(configOpts):
    from Error import raiseError
    # Read a dust model from an options dictionary.
    # Depending on the dust model, certain keywords are expected to be present.
    # If it's a map then we read it and return it
    # Additional dust models can be implemented by adding code here and in
    # determine_dust
    implemented = ['const', 'slope', 'map']
    model_type = configOpts['dusttype']
    if model_type not in implemented:
        raiseError('Fatal error: no valid dust model specified.')
    if model_type == 'const':
        return model_type, [configOpts['dustval'], configOpts['dusterr']]
    elif model_type == 'slope':
        params = [configOpts['slope'], configOpts[
            'offset'], configOpts['slopeerror']]
        return model_type, params
    # return the filename of the dust map and solar metallicity for conv.
    elif model_type == 'map':
        return model_type, [open_image(configOpts['mapname'])[0], open_image(configOpts['mapname'])[1],
                            configOpts['solar'], open_image(configOpts['maperror'])[1]]
        # The error map is assumed to have the same WCS header as the
        # metallicity map


def read_array(filename, dtype, separator=',', comments='#'):
    # adapted from the scipy cookbook (to ignore comment lines)
    # Read a file with an arbitrary number of columns.
    # The type of data in each column is arbitrary
    # It will be cast to the given dtype at runtime
    data = [[] for dummy in xrange(len(dtype))]
    for line in open(filename, 'r'):
        if line[0] != comments:
            fields = line.strip().split(separator)
            for i, number in enumerate(fields):
                data[i].append(number)
    # create an array of numpy arrays or rec.array won't work
    for i in xrange(len(dtype)):
        data[i] = np.cast[dtype[i]](data[i])
                # print data
                # print np.rec.array(data, dtype=dtype)
    return np.rec.array(data, dtype=dtype)


def read_coords(filename, separator=',',  comments='#'):
    # reads a list of coordinates
    columns = np.dtype([('PDRID', 'int32'), ('RA', '|S14'), ('DEC', '|S14')])
    data = read_array(filename, columns, separator, comments)
                # print data
                # print np.rec.array(data, dtype=columns)
    return np.rec.array(data, dtype=columns)


def open_image(filename):
    # opens a fits file and returns an astWCS.WCS object and a FITSimage object that also has a header etc. (Prev. version see pdrLib.py.4)
    # reverted to original version since we don't use Kapteyn anymore
    # the numpy array should be image (raw image)
    #image = maputils.FITSimage(filename)
    # return astWCS.WCS(image.hdr, mode="pyfits"), image
    image = pf.open(filename)
    header = image[0].header

    keys = ['CTYPE', 'CDELT', 'CRVAL', 'CRPIX', 'CROTA']

    if header['NAXIS'] != 2:
        try:
            # print 'Attempting to extract image' # Replace with proper log!
            # print 'If this fails, please pre-flatten the image to 2D'
            slicestring = ''
            for n in range(3, header['NAXIS'] + 1):
                #Remove all keywords related to axes > 2
                slicestring += '0,'
                for k in keys:
                    if k + str(n) in header:
                        del header[k + str(n)]
            header['NAXIS'] = 2
            #Select only the RA and DEC axes
            data = eval('image[0].data['+slicestring+':,:]')
            #Remember that pyfits axes are read inverted; e.g. STOKES, MOM0, DEC, RA
            image[0].header = header
            image[0].data = data
        except:
            raiseError('Trying to flatten image %s to 2-D failed. ' % filename +
                       'Please, do it manually.')
    else:
        data = image[0].data

    return pw.WCS(header), data, image


def logdump(filename, line):
    with open(filename, 'a') as f:
        f.write(line)

""" Various PDR method calculations """


def PDRmodel(flux, ext, D_gal, dd0, rho, NHI):
    # We now expect Dgal in pc, and rho in pc!
                # Nowhere do we output G0, whereas people want to know. Quick hack to print it here, might return it as a tuple
    # scale NHI; and calculate a number-loss friendly conversion factor

    N = NHI * 1e-21 / 0.78
    # As in Heiner et al. 2008a
    conv = 106 * 0.85 * 10 ** (ext / 2.5) / 0.264 * (D_gal / 1e5) ** 2
    # print "  conv ", conv
    # This is duplicate now that we calculate G0 earlier in the pipeline!
    G0 = 0.85 * flux * 10 ** (ext / 2.5) * (D_gal / rho) ** 2. / 2.64e-6
    # print "G0: ", G0 #to check if the calculated G0s are identical
    #"Original" model:
    # return (conv * flux*1e15) / (rho**2 * m.sqrt(dd0) * (m.exp(N*dd0) - 1)) #n
    # Model with Heaton's improvement (Heaton 2009, Eq. 12):
    # avoid division by zero
    # rewrote expression because denominator can get very large for strange dd0
    denominator = np.log10(rho ** 2 * m.sqrt(dd0) * (m.exp(N * dd0) - 1))
    if (denominator == 0.):
        return 0.
    else:
        numerator = np.log10(conv * flux * 1e15 * (dd0 ** 0.2))
        # return (conv * flux*1e15) * (dd0**0.2) / denominator #n
        return 10 ** (numerator - denominator)


def Gnaught(flux, D_gal, rho, ext):  # D_gal must be in pc
    return 0.85 * flux * 10 ** (ext / 2.5) * (D_gal / rho) ** 2. / 2.64e-6


def contrast(flux_raw, mean_at_r, rho, D_gal):
    #!! mean_at_r is already per arcsec^2 (see generic_flux)
    # arc_to_sterad in arcsec^2 / sterad; integration over 'pi' sterad
    arc_to_sterad = 42545170300
    fuv_at_HI = flux_raw * (D_gal / rho) ** 2  # D_gal must be in pc
    bg_level = mean_at_r * arc_to_sterad * m.pi
    # handle cases where mean_at_r = 0 by setting a low contrast (likely to
    # get filtered)
    if (bg_level == 0):
        return 1e-6
    else:
        return fuv_at_HI / bg_level  # ratio source/background


def errors(n, flux, sF, rho, srho, dd0, sdd0, NHI, sNHI, ext, logfile=False):
    # Formulas described in Heiner et al. 2008a; updated Heaton09
    # All errors come in as absolute errors
    #n = PDRmodel(F, ext, dd0, rho, NHI)
    # Some reworking to minimize number loss and calculation
    N = NHI / 0.78e21
    sN = sNHI / 0.78e21
    if (dd0 == 0. or N == 0.):  # avoid division by zero; set expfrac high
        expfrac = 1e6
    else:
        expfrac = m.exp(N * dd0) / (m.exp(N * dd0) - 1)
    if (flux != 0.):
        one = (sF / flux) ** 2
    else:
        # avoid division by zero and artificially inflate the error in such a
        # case
        one = 1e6
    two = 4 * (srho / rho) ** 2
    three = (-0.3 / dd0 - N * expfrac) ** 2 * sdd0 ** 2
    four = (-dd0 * expfrac * sN) ** 2
    sig_n = n * m.sqrt(one + two + three + four)
    if (n == 0.):
        rel_err = 0.
    else:
        rel_err = sig_n / n * 100
    if (logfile):
        logdump(logfile, "{} {} {} {} {} {} {}\n".format(
            one, two, three, four, n, sig_n, rel_err))
    # JSH 2012-06-04 with strange values of dd0 can end up with nan
    # which would be appropriate for a value of which we cannot determine
    # error
    return sig_n


def determine_dust(scenario, p):
    # Various scenarios to determine d/d0 are recognized through parameter array 'p'
    # Uses separation(...)
    # JSH 2010-02-06; 2010-11-24 (removed calculation of R)
    # JSH 2012-06-04: limit 12+log(O/H) to 9.5 in map due to extreme values (n < 0.014)
    # Issue: when passing a filename slope still gets parsed and returns an error.  << no it was a calling error
    # Workaround: use 'try' and 'except' in slope so it doesn't get parsed before it is needed
    # Wishlist: replace 'p' by more transparent solution

    def const():
        # Always return the same dd0; the array p is standard but p[2] will
        # contain the constant
        return p[2][0], p[2][1]  # value, relative error

    def slope():
        # dd0(R) = 10^(a*R + C); a and C are generally negative
        # p: [[RA,DEC],[center_RA,center_DEC], [a, C, err], [pa, incl, D_gal(kpc), Rgal], astWCS header]
        # Determine galactocentric radius here from coord & center
        #R = separation(p[0][0], p[0][1], p[1][0], p[1][1], p[3][0], p[3][1], p[3][2], p[4])
        # print p[0][0], p[0][1], p[1][0], p[1][1], R
        #dd0 = 0
        # try:
        # print(float(p[2][0]), p[3][3], float(p[2][1]))
        # print 10 ** (float(p[2][0]) * p[3][3] + float(p[2][1]))
        dd0 = 10 ** (float(p[2][0]) * p[3][3] + float(p[2][1]))
        # constant relative error read from inputConfig.dat
        sdd0 = p[2][2] * dd0
        # It is more convenient to estimate a relative error in this case, but the code expects
        # an absolute error
        # except:
        #	print "Error in dust model slope."
        return dd0, sdd0

    # def anyfunction(): #easy to implement a custom one with its own
    # parameter array

    def map():
        # p: [[RA,DEC],[center_RA,center_DEC], [header, image, solar, err_img],
        # [pa, incl, D_gal(kpc), Rgal], astWCS header, logger]
        # 2013-11-13: updated p to include logger

        # Returns the pixel value at a given RA, DEC
        # p: [[RA, DEC], astWCS header, numpy image array]
        # print "Converting to d/d0 with solar 12+log(O/H)={0}".format(p[2][2])

        # 2013-11-13 re-write to use astropy only
        #            todo: deal with logger (needs to be passed)
        #                  perhaps as last element of p
        from astropy import coordinates as coord
        from astropy import wcs

        logger = p[5]

        #Ultimately we should update the use of p so we don't have to
        #convert locally
        coordstring = p[0][0].tostring() + ' ' + p[0][1].tostring()
        logger.write(coordstring, newLine=False)
        location = coord.ICRSCoordinates(coordstring)

        w = p[4] #contains the header object
        r, d = (location.ra.degrees, location.dec.degrees)
        x, y = np.round(w.wcs_world2pix([[r, d]], 0))[0] #second argument is base 0; set to 1 if you want DS9-like values
        ohdata = p[2][1]
        sohdata= p[2][3]
        ymax, xmax = ohdata.shape
        #Need to check that the coordinates are in fact inside the
        #image since Python will accept negative indices
        oh = 0
        soh = 0.0001 #bit of a hack to distinguish between hardcoded soh
        #If an oh value is found soh will also be non-zero
        if (x < 0 or x >= xmax) or (y < 0 or y >= ymax):
            logger.write("Warning: coordinate outside map; setting 0", newLine=False)
        else:
            oh = ohdata[y, x]
            if np.isnan(oh) or oh == 0:
                try:
                    # If there is no value at this location we will take the
                    # highest value within a small box and warn about it
                    # (ignore NaNs)
                    oh = np.nanmax(p[2][1][y - 2:y + 3, x - 2:x + 3])
                except:
                    logger.write("Tried to scan outside image; setting oh to 0", newLine=False)
                if (m.isnan(oh) or oh == 0):
                    logger.write("Warning: no value found", newLine=False)
                else:
                    soh = p[2][3][y - 2:y + 3, x - 2:x + 3].ravel()[
                        np.nanargmax(p[2][1][y - 2:y + 3, x - 2:x + 3])]
                    # Note that the proper index applies to a subselection of
                    # both images ^
                    logger.write("Read 12+log(O/H) = {} +/- {} -- maximum value nearby".format(oh, soh), newLine=False)
            else:
                soh = p[2][3][y, x]
                logger.write("Read 12+log(O/H) = {} +/- {}".format(oh, soh), newLine=False)

        if (oh > 9.5):
            logger.write("Warning: OH > 9.5 rejected, set to 0", newLine=False)
            oh = 0

        dd0 = 10 ** (oh - p[2][2])
        # absolute error; all errors absolute to avoid confusion
        sdd0 = dd0 * p[2][2] * soh / 10.
        # sdd0 = p[2][2] * soh / 10. #relative error
        # a value of oh = 0 will result in artificially high n
        return dd0, sdd0

    # print scenario, p

    options = {'const': const,
               'slope': slope,
               'map': map}

    return options[scenario]()

""" Auxilliary functions """


def records(table, labels):  # generic numpy.rec.array wrapper to label data
        # expect a table of the form [[], [], ...]
        # labels as 'name1, name2, ...'
    for i in range(len(table)):
        table[i] = np.array(table[i])

    return np.rec.array(table, names=labels)


""" Example presentation routines """


def plot_n(R, logn):
        # Custom ranges to filter out values made artificially high
    plt.plot(R, logn, 'ko', markersize=7, markeredgewidth=2)
    plt.xlabel(r"$\rm{R/R_{25}}$", fontsize=17)
    plt.ylabel(r"$\rm{log\ (n\ cm^{-3})}$", fontsize=17)
    plt.axis([0, 1, -3, 3])
    plt.show()
    return 0


def plot_all(X, Y, Xlabel, Ylabels, loglist):

    y_offset = 0.02  # 0.01

    no_plots = np.shape(Y)[0]
    for i in range(no_plots):
        panel = plt.subplot(no_plots, 1, i + 1)
        position = panel.get_position()
        pos = position.get_points()
        #pos[0][0] = 0.13
        pos[0][1] = pos[0][1] + y_offset
        position.set_points(pos)
        panel.set_position(position)

        if loglist[i]:
            plt.semilogy(X, Y[i], 'ko')
        else:
            plt.plot(X, Y[i], 'ko')
        plt.ylabel(Ylabels[i])
    plt.xlabel(Xlabel)
    plt.show()
    return 0


def ntot_ds9(filename, data):
    regfile = open(filename, 'w')
    regfile.write(
        "#Resulting total hydrogen volume densities per patch\n#Radius = 2*(floor(log10(ntot))+3)\n")
    regfile.write("global color = green\nfk5\n")
    for i in range(np.size(data.ntot)):
        # Values lower than 0.01 get same radius
        radius = 2 * (np.floor(np.log10(data.ntot[i])) + 3)
        if (radius < 1):
            radius = 1.
        ra = data.pRA[i].strip().replace(' ', ':')
        dec = data.pDEC[i].strip().replace(' ', ':')
        regfile.write("circle({0},{1},{2}\")\n".format(ra, dec, radius))
    regfile.close()
    return 0
