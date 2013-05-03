#!/usr/bin/env python
# encoding: utf-8

#Collection of functions for various parts of the PDR method 'pipeline' script pipeline.py
#Depends on astLib (0.4.0), Pyfits, NumPy & SciPy, (kapteyn), matplotlib, sextutils
#JSH 2010-02-01
#    2011-07-14
#    2012-02-14~latest edit (needed to fix a bug on Fedora with formatting)
#    2012-02-15 edits to avoid divisions by zero (should disqualify the measurement)
#                 also, fix errors when going outside the dust map
# & José Sanchez-Gallego
#    2012-06-04 math range errors for high values of 12+log(O/H)

# Fedora bug: https://bugzilla.redhat.com/show_bug.cgi?format=multiple&id=654269
#   Workaround is to 'dumb down' the formatting statements and do a type cast in advance
#   It seems to happen with 'x' and 'd' formatters
#   Another workaround is to replace the 'd' with 'g', which has the desired effect

#pdrLib is now organized to contain José's routines as well, seperate from this file.

#Wishlist: classes like coordinate lists with built in options to query hms dms, or decimal etc.
#		because now we have to do all kinds of conversions first
#		can get rid of passing the header explicitely all the time

#from __future__ import with_statement #Needed for logdump's 'with' use in Python 2.5 only?
import pdb, sys #debugger
import numpy as np
from impPyFits import *
from string import strip
from astLib import astCoords, astWCS #astlib.sourceforge.net coordinate conversions & calculations
from matplotlib import pylab as plt #for flux specifically and plot routines in general
from os import system #to call SExtractor
from sextutils import se_catalog #Ferguson's code to read a SExtractor catalog

""" Default file operations """


def read_dust(configOpts):
  #Read a dust model from an options dictionary.
  #Depending on the dust model, certain keywords are expected to be present.
  #If it's a map then we read it and return it
  #Additional dust models can be implemented by adding code here and in determine_dust
  implemented = ['const', 'slope', 'map']
  model_type = configOpts['dusttype']
  if model_type not in implemented:
    print "Fatal error: no valid dust model specified."
    exit(0)
  if model_type == 'const':
    return model_type, [configOpts['dustval'], configOpts['dusterr']]
  elif model_type == 'slope':
    params = [configOpts['slope'], configOpts['offset'], configOpts['slopeerror']]
    return model_type, params
  elif model_type == 'map': #return the filename of the dust map and solar metallicity for conv.
    return model_type, [open_image(configOpts['mapname'])[0], open_image(configOpts['mapname'])[1], configOpts['solar'], open_image(configOpts['maperror'])[1]]
    #The error map is assumed to have the same WCS header as the metallicity map

def read_array(filename, dtype, separator=',', comments='#'):
  #adapted from the scipy cookbook (to ignore comment lines)
  #Read a file with an arbitrary number of columns.
  #The type of data in each column is arbitrary
  #It will be cast to the given dtype at runtime
  data = [[] for dummy in xrange(len(dtype))]
  for line in open(filename, 'r'):
    if line[0] != comments:
      fields = line.strip().split(separator)
      for i, number in enumerate(fields):
        data[i].append(number)
  for i in xrange(len(dtype)): #create an array of numpy arrays or rec.array won't work
    data[i] = np.cast[dtype[i]](data[i])
        #print data
        #print np.rec.array(data, dtype=dtype)
  return np.rec.array(data, dtype=dtype)

def read_coords(filename, separator=',',  comments='#'):
  #reads a list of coordinates
  columns = np.dtype([('PDRID', 'int32'), ('RA', '|S14'), ('DEC', '|S14')])
  data = read_array(filename, columns, separator, comments)
        #print data
        #print np.rec.array(data, dtype=columns)
  return np.rec.array(data, dtype=columns)

def open_image(filename):
  #opens a fits file and returns an astWCS.WCS object and a FITSimage object that also has a header etc. (Prev. version see pdrLib.py.4)
  #reverted to original version since we don't use Kapteyn anymore
  #the numpy array should be image (raw image)
  #image = maputils.FITSimage(filename)
  #return astWCS.WCS(image.hdr, mode="pyfits"), image
  image = pf.open(filename)
  return astWCS.WCS(image[0].header, mode="pyfits"), image[0].data

def logdump(filename, line):
  with open(filename, 'a') as f:
    f.write(line)

""" Various PDR method calculations """

def PDRmodel(flux, ext, D_gal, dd0, rho, NHI):
  #We now expect Dgal in pc, and rho in pc!
        #Nowhere do we output G0, whereas people want to know. Quick hack to print it here, might return it as a tuple
  #scale NHI; and calculate a number-loss friendly conversion factor
  N = NHI * 1e-21 / 0.78
  conv = 106 * 0.85 * 10**(ext/2.5) / 0.264 * (D_gal/1e5)**2 #As in Heiner et al. 2008a
  #print "  conv ", conv
  #This is duplicate now that we calculate G0 earlier in the pipeline!
  G0 = 0.85 * flux * 10**(ext/2.5) * (D_gal/rho)**2. / 2.64e-6
  #print "G0: ", G0 #to check if the calculated G0s are identical
  #"Original" model:
  #return (conv * flux*1e15) / (rho**2 * m.sqrt(dd0) * (m.exp(N*dd0) - 1)) #n
  #Model with Heaton's improvement (Heaton 2009, Eq. 12):
  #avoid division by zero
  #rewrote expression because denominator can get very large for strange dd0
  denominator = np.log10(rho**2 * m.sqrt(dd0) * (m.exp(N*dd0) - 1))
  if (denominator == 0.):
    return 0.
  else:
    numerator = np.log10(conv * flux*1e15 * (dd0**0.2))
    #return (conv * flux*1e15) * (dd0**0.2) / denominator #n
    return 10**(numerator - denominator)


def Gnaught(flux, D_gal, rho, ext): #D_gal must be in pc
  return 0.85 * flux * 10**(ext/2.5) * (D_gal/rho)**2. / 2.64e-6

def e_angle_rad(x1, y1, x2, y2, pa):
  	#Returns angles between two vectors, taking pa into account
	#Since the output is used in a (symmetrical) ellipse equation,
	#we don't care about the quadrants
	#Assume that the first set is the center, but does it matter?
	#need to verify, but currently first set is always center

	#1. 'center' unit vector, rotated by pa
	xc = m.cos(m.radians(pa+90))
	yc = m.sin(m.radians(pa+90))
	#2. construct vector from first coord to second coord
	xp = x2 - x1
	yp = y2 - y1
	#3. calculate angle between 1. and 2.
	lc = m.sqrt(xc**2. + yc**2.)
	lp = m.sqrt(xp**2. + yp**2.)
	#print x1, y1, x2, y2, xc, yc, xp, yp, "input values\n", lc, lp, "lc lp"
	x = (xc*xp+yc*yp)/(lc*lp)
	y = m.sqrt(1-x**2.)
	#print x, y, m.atan2(y,x), m.degrees(m.atan2(y,x)), "angle between vectors"
	return(m.atan2(y,x));

def separation(ra1, dec1, ra2, dec2, pa, incl, D_gal, wcsheader):
  logfile = "rhodump.log"
  logdump(logfile, "{0} {1} {2} {3}\n".format(ra1, dec1, ra2, dec2))

  #Assuming ra & dec come in as strings as "hh mm ss.sss"
  #1. calculate separation angle between the two coordinates, giving us r_raw
  #Strip: can't have leading or trailing whitespaces
  r1, d1 = astCoords.hms2decimal(strip(ra1), " "), astCoords.dms2decimal(strip(dec1), " ")
  r2, d2 = astCoords.hms2decimal(strip(ra2), " "), astCoords.dms2decimal(strip(dec2), " ")
  angle = astCoords.calcAngSepDeg(r1, d1, r2, d2)
  r_raw = D_gal * m.radians(angle)
  #print r_raw, "r_raw"

  #2. correct for pa, i by deprojecting radius
  #For this we need x and y coordinates through the wcs header
  x1, y1 = wcsheader.wcs2pix(r1, d1)
  x2, y2 = wcsheader.wcs2pix(r2, d2)
  e_a = e_angle_rad(x1, y1, x2, y2, pa) #get the angle
  esqd = 1 - (m.cos(m.radians(incl)))**2. #deproject
  sqrtf = m.sqrt((1-esqd)/(1-esqd*m.cos(e_a)**2.))
  #print esqd, m.cos(e_a), sqrtf, "esqd cos(angle) sqrtf"
  #calculate separation in parsec from angle and ellipse deprojection

  logdump(logfile, "{0} {1} {2}\n".format(r_raw, sqrtf, r_raw/sqrtf))

  #Append to a ds9-style file - replace our ' ' with ':'
  r1s = ra1.strip().replace(' ', ':')
  d1s = dec1.strip().replace(' ', ':')
  r2s = ra2.strip().replace(' ', ':')
  d2s = dec2.strip().replace(' ', ':')
  logdump("rholog.reg", "circle({0},{1},{2}\")\n".format(r1s, d1s, r_raw/sqrtf/D_gal*3600))
  logdump("rholog.reg", "point({0},{1})\n".format(r2s,d2s))

  return r_raw / sqrtf #in units of D_gal (e.g. kpc)

def contrast(flux_raw, mean_at_r, rho, D_gal):
  #!! mean_at_r is already per arcsec^2 (see generic_flux)
  #arc_to_sterad in arcsec^2 / sterad; integration over 'pi' sterad
  arc_to_sterad = 42545170300
  fuv_at_HI = flux_raw*(D_gal/rho)**2   #D_gal must be in pc
  bg_level = mean_at_r * arc_to_sterad * m.pi
  #handle cases where mean_at_r = 0 by setting a low contrast (likely to get filtered)
  if (bg_level == 0):
    return 1e-6
  else:
    return fuv_at_HI / bg_level #ratio source/background

def errors(n, flux, sF, rho, srho, dd0, sdd0, NHI, sNHI, ext, logfile = False):
  #Formulas described in Heiner et al. 2008a
  #All errors come in as absolute errors
  #n = PDRmodel(F, ext, dd0, rho, NHI)
  #Some reworking to minimize number loss and calculation
  N = NHI / 0.78e21
  sN = sNHI / 0.78e21
  if (dd0 == 0. or N == 0.): #avoid division by zero; set expfrac high
    expfrac = 1e6
  else:
    expfrac = m.exp(N*dd0) / (m.exp(N*dd0) - 1)
  if (flux != 0.):
    one = (sF / flux)**2
  else:
    one = 1e6 #avoid division by zero and artificially inflate the error in such a case
  two = 4 * (srho / rho)**2
  three = (-1/(2*dd0) - N*expfrac)**2 * sdd0**2
  four = (-dd0 * expfrac * sN)**2
  sig_n = n * m.sqrt(one+two+three+four)
  if (n == 0.):
    rel_err = 0.
  else:
    rel_err = sig_n / n * 100
  if (logfile):
    logdump(logfile, "{} {} {} {} {} {} {}\n".format(one, two, three, four, n, sig_n, rel_err))
  #JSH 2012-06-04 with strange values of dd0 can end up with nan
  #  which would be appropriate for a value of which we cannot determine error
  return sig_n

def determine_dust(scenario,p):
  #Various scenarios to determine d/d0 are recognized through parameter array 'p'
  #Uses separation(...)
  #JSH 2010-02-06; 2010-11-24 (removed calculation of R)
  #JSH 2012-06-04: limit 12+log(O/H) to 9.5 in map due to extreme values (n < 0.014)
  #Issue: when passing a filename slope still gets parsed and returns an error.  << no it was a calling error
  #Workaround: use 'try' and 'except' in slope so it doesn't get parsed before it is needed
  #Wishlist: replace 'p' by more transparent solution

  def const():
    #Always return the same dd0; the array p is standard but p[2] will contain the constant
    return p[2][0], p[2][1] #value, relative error

  def slope():
    #dd0(R) = 10^(a*R + C); a and C are generally negative
    #p: [[RA,DEC],[center_RA,center_DEC], [a, C, err], [pa, incl, D_gal(kpc), Rgal], astWCS header]
    #Determine galactocentric radius here from coord & center
    #R = separation(p[0][0], p[0][1], p[1][0], p[1][1], p[3][0], p[3][1], p[3][2], p[4])
    #print p[0][0], p[0][1], p[1][0], p[1][1], R
    #dd0 = 0
    #try:
    dd0 = 10**(float(p[2][0])*p[3][3] + float(p[2][1]))
    sdd0 = p[2][2] * dd0 #constant relative error read from inputConfig.dat
    #It is more convenient to estimate a relative error in this case, but the code expects
    #an absolute error
    #except:
    #	print "Error in dust model slope."
    return dd0, sdd0

  #def anyfunction(): #easy to implement a custom one with its own parameter array

  def map():
    #p: [[RA,DEC],[center_RA,center_DEC], [header, image, solar, err_img], [pa, incl, D_gal(kpc), Rgal], astWCS header]

    #Returns the pixel value at a given RA, DEC
    #p: [[RA, DEC], astWCS header, numpy image array]
    #print "Converting to d/d0 with solar 12+log(O/H)={0}".format(p[2][2])

    r, d = astCoords.hms2decimal(strip(p[0][0]), " "), astCoords.dms2decimal(strip(p[0][1]), " ")
    oh = soh = 0
    if p[2][0].coordsAreInImage(r, d):
      x, y = np.round(p[2][0].wcs2pix(r, d))
      #astLib returns indices for the image array; 1 lower in value than DS9 x, y
    #Assume the image is the first entry in the hdulist
    #PyFits Handbook p. 16: slow axis first, last axis last,

    #!At this time we assume a 12+log(O/H) map so the conversion is hardcoded and should be parameters!!:
      oh =  p[2][1][y,x]
      if (m.isnan(oh) or oh==0):
        #Even if p[2][0] is in image, this box might not be
        try:
          #If there is no value at this location we will take the highest value within a small box and warn about it (ignore NaNs)
          oh = np.nanmax(p[2][1][y-2:y+3,x-2:x+3])
        except:
          print "Tried to scan outside image; setting oh to 0"
        if (m.isnan(oh) or oh==0):
          print "Warning: no value found"
        else:
          soh = p[2][3][y-2:y+3,x-2:x+3].ravel()[np.nanargmax(p[2][1][y-2:y+3,x-2:x+3])]
          #Note that the proper index applies to a subselection of both images ^
          print "Read 12+log(O/H) = {} +/- {} -- maximum value nearby".format(oh, soh)
      else:
        soh = p[2][3][y,x]
        print "Read 12+log(O/H) = {} +/- {}".format(oh, soh)
    else: #coordinate is not on the map
        print "Warning: coordinate outside map; setting 0"
    #R = separation(p[0][0], p[0][1], p[1][0], p[1][1], p[3][0], p[3][1], p[3][2], p[4])

    if (oh > 9.5):
      print "Warning: OH > 9.5 rejected, set to 0"
      oh = 0

    dd0 = 10**(oh - p[2][2])
    sdd0 = dd0 * p[2][2] * soh / 10. #absolute error; all errors absolute to avoid confusion
    #sdd0 = p[2][2] * soh / 10. #relative error
    #a value of oh = 0 will result in artificially high n
    return dd0, sdd0

  #print scenario, p

  options = {'const': const,
             'slope': slope,
             'map': map}

  return options[scenario]()

""" Auxilliary functions """
def records(table, labels): #generic numpy.rec.array wrapper to label data
    #expect a table of the form [[], [], ...]
    #labels as 'name1, name2, ...'
    for i in range(len(table)):
        table[i] = np.array(table[i])

    return np.rec.array(table, names=labels)

def generic_flux(filename, coordsPix, maxrad_pix, pix_size = 1.0, contrast = 0.5, fluxerror = 0.3, logfile = False):

    from calcFlux import calcFluxMulti

    genericFluxResults = [] # The list with the final results to return

    # Calls calcFluxMulti an gets the areas and fluxes for each position and aperture
    # from 1 to maxrad_pix
    fluxes = calcFluxMulti(filename, coordsPix, 1, maxrad_pix, step=1)

    for pos in range(coordsPix.shape[0]): # iterates over the different positions

        fluxlist = []
        avglist = []
        mean_at_r = 0; aperture = 0; bgflux = 0; cumulflux = 0; netflux = 0; sflux = 0

        for n in range(maxrad_pix): #iterate over rings of width 1 pixel

            r = n + 1
            fluxlist.append([fluxes[pos][5][n], fluxes[pos][6][n]])
            # print fluxlist
            if n > 0:
                #incremental flux divided by ring surface area
                ring_avg = (fluxlist[n][1] - fluxlist[n-1][1]) / (fluxlist[n][0] - fluxlist[n-1][0])
                avglist.append(ring_avg) #current avglist index will be n-1 since it starts at n=1

                if (logfile):
                    logdump(logfile, "{0:3g} {1:7.2e} {2:8.3e} {3:7.5e}\n".format(r, fluxlist[n][0], fluxlist[n][1], ring_avg))

            if (n > 1) and (not mean_at_r):
                #Check if the mean value has dropped to at least 50% (contrast) of the initial value (arbitrary)
                #  note that if the actual max value is much higher, 50% is not indicative of any kind of source contrast
                #Then check if a local minimum was encountered.
                #This could create problems with very poor source contrast, but those may need to be filtered out earlier on
                #For example: need to drop to at least half the initial mean value first
                if (avglist[n-1] < (avglist[0] * contrast)):
                    if (avglist[n-1] > avglist[n-2]): #list is 1 shorter than n
                        mean_at_r = avglist[n-2] #per pixel
                        aperture = (r-1) * pix_size #in arcsec and assuming square pixels
                        bgflux = mean_at_r * fluxlist[n-2][0] #flux/pixel*area[in pixels]
                        cumulflux = fluxlist[n-2][1]
                        netflux = cumulflux - bgflux
                        mean_at_r = mean_at_r / (pix_size*pix_size) #!need to report per arcsec^2
                        sflux = fluxerror*netflux #flux error as in config file (fixed relative error)

                        if (logfile):
                          logdump(logfile, "#First local minimum registered\n")
                        #print aperture, r-1

        if (not mean_at_r): #we exit leaving all the values at zero
            print "Warning: minimum not found! (pos {0})".format(pos)
            if (logfile):
                logdump(logfile, "#Warning: minimum not found!\n")

        genericFluxResults.append([mean_at_r, aperture, bgflux, cumulflux, netflux, sflux])

    return genericFluxResults


# fuv_flux
def fuv_flux(imageName, header, coords, configOpts, verbose=False):
  #Wishlist: use qphot to get all the annuli in one go
  #configOpts used: FUVscale, contrast, maxrad, uv_log
  #Determines the (FUV) flux at a specified coordinate
  #A net flux is calculated based on the a local minimum in ring-averaged flux values per pixel
  #The maximum radius can't be bigger than half the box size (not checked)
  #Need to implement: option to accept first local minimum or go to next

  logfile = configOpts['uv_log']

  columns = np.dtype([('PDRID', 'int32'), ('RA', '|S14'), ('DEC', '|S14'), ('aperture', 'float'), ('mean_at_r', 'float'), ('bgflux', 'float'), ('cumulflux', 'float'), ('netflux', 'float'), ('sflux', 'float')])
  results = []#[[] for dummy in xrange(len(columns))]
  pix_size = [header.getXPixelSizeDeg() * 3600, header.getYPixelSizeDeg() * 3600]
  #pix_area = pix_size[0] * pix_size[1]
  logdump(logfile, "#Number of coordinates: {0}\n".format(len(coords)))
  logdump(logfile, "#R(pix) Area (pix) Flux (image units) Flux/pix\n")

  #Determine flux for all coordinates

  maxrad_pix = int(round(configOpts['maxrad']/pix_size[0])) # Global maximum radius (in pixels)
  coordsPix = np.zeros([len(coords), 2], np.float)

  for p in range(len(coords)):
    #Define circle center and check pixel size
    x, y = np.round(header.wcs2pix(astCoords.hms2decimal(strip(coords[p].RA), " "), astCoords.dms2decimal(strip(coords[p].DEC), " "))) #x, y not integers so int(x) will trunc not round
    #x, y should be rounded values or requesting dummy[y,x] will lead to truncation errors
    logdump(logfile, "#Source {0} at pixel position: {1}, {2}\n".format(coords[p].PDRID, x, y))
    coordsPix[p, :] = [x, y]

  genericFluxResults =  generic_flux(imageName, coordsPix, maxrad_pix, pix_size[0],
    configOpts['contrast'], configOpts['fluxerror'], logfile = logfile)

  for p in range(len(coords)):

    mean_at_r, aperture, bgflux, cumulflux, netflux, sflux = genericFluxResults[p]

    #Scale the results -- for generic use set FUVscale to 0
    scale = configOpts['FUVscale'] # Converts from cps to flux units
    mean_at_r *= scale
    bgflux *= scale
    cumulflux *= scale
    netflux *= scale
    sflux *= scale

    #Store the results
    if verbose:
      print "local minimum at", mean_at_r, aperture
      print "bgflux", bgflux
      print "cumul flux", cumulflux
      print "net flux", netflux

    results.append([coords[p].PDRID, coords[p].RA, coords[p].DEC, aperture, mean_at_r, bgflux, cumulflux, netflux, sflux])

  return np.rec.array(results, dtype=columns)

### HI patches ###

# Run SExtractor, extract the results
# SExtractor requires files default.conv default.nnw pdrSExout.param and HIclumps.sex
# pdrLib requires only that pdrSEx.param include ALPHAPEAK_J2000 and DELTAPEAK_J2000

#Call SExtractor - this procedure is simply a specialized wrapper to SExtractor
#Need SExtractor in path
#It needs to be called for every postage stamp HI region
def call_SEx(fitsfile):
  #The HIclumps.sex config file comes with this software. It does not need to be customized further.
  #Output catalog: HIclumps.cat
  err = system("sex " + fitsfile + " -c HIclumps.sex")
  if (err):
    print "Error calling SExtractor."
    return 1
  #If no clumps were found, the .cat file will only contain commented lines
  return "HIclumps.cat" #This file only needs to be hard-coded here and in HIclumps.sex

#Call Fergusons's library to read a SExtractor catalog,
#produced previously to select HI patches in an individual candidate PDRs
def read_secat(catfile, fitsfile, logfile=False, wcs=False):
  patchtable = [[] for dummy in xrange(3)]  #ra, dec, NHI
  #Must catch empty catalog file:
  try: #Call will fail if the catalog file is empty or does not have the expected format
    catalog = se_catalog(catfile, readfile=True, preserve_case=False)
  except:
    print "Warning: no clumps found. Creating empty entry."
    patchtable[0].append(0), patchtable[1].append(0)
    patchtable[2].append(0)
    return patchtable

  records = len(catalog.alphapeak_j2000)
  header, image = open_image(fitsfile)
  for n in range(records):
    a, d = catalog.alphapeak_j2000[n], catalog.deltapeak_j2000[n]
    x, y = header.wcs2pix(a, d)
    x = int(np.round(x))
    y = int(np.round(y))
    maxHI = image[y, x]
    print "a,d: ", a, d, "x,y: ", x, y, "value: ", maxHI #no offset? Verify!!
    patchtable[0].append(astCoords.decimal2hms(a, " "))
    patchtable[1].append(astCoords.decimal2dms(d, " "))
    patchtable[2].append(maxHI)
    #patchtable[3].append(0.1)    #error now calculated after calling this routine
    #print "Warning: hardcoding sNHI at 0.1"
    logdump(logfile, "{} {} {}\n".format(a, d, maxHI))

  return patchtable #calling function administrates PDRID


""" Example presentation routines """

def plot_n(R, logn):
	#Custom ranges to filter out values made artificially high
  plt.plot(R, logn, 'ko', markersize=7, markeredgewidth=2)
  plt.xlabel(r"$\rm{R/R_{25}}$", fontsize = 17)
  plt.ylabel(r"$\rm{log\ (n\ cm^{-3})}$", fontsize = 17)
  plt.axis([0, 1, -3, 3])
  plt.show()
  return 0

def plot_all(X, Y, Xlabel, Ylabels, loglist):

  y_offset = 0.02 #0.01

  no_plots = np.shape(Y)[0]
  for i in range(no_plots):
    panel = plt.subplot(no_plots, 1, i+1)
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
  regfile.write("#Resulting total hydrogen volume densities per patch\n#Radius = 2*(floor(log10(ntot))+3)\n")
  regfile.write("global color = green\nfk5\n")
  for i in range(np.size(data.ntot)):
    radius = 2*(np.floor(np.log10(data.ntot[i])) + 3) #Values lower than 0.01 get same radius
    if (radius < 1):
      radius = 1.
    ra = data.pRA[i].strip().replace(' ', ':')
    dec = data.pDEC[i].strip().replace(' ', ':')
    regfile.write("circle({0},{1},{2}\")\n".format(ra, dec, radius))
  regfile.close()
  return 0

