#! /usr/bin/env python
# encoding: utf-8
#
# Usage: PMAP.py [configfile]
# This is the main pipeline script aiming to provide a one-stop PDR method solution.
# Several steps will be skipped if certain files already exist from previous runs.
# A quick plot is produced with the results.
# 2011-08-04: now uses python2.7 because of Pyraf (for flux calculation)
#
# JSH 2010-10-28
#
#     previous edit: 2011-10-26
#     last edit: 2012-02-14 -- fedora bug w/ formatting see pdrLib.py also
#
#     New edit: 2013-04-13 by José R. Sánchez-Gallego
#         -- General clean-up and optimisation
#
#     2013-04-26 update logging to use fancyPrint << not tested

# import pdb  # debugger
import pdrLib as pdr
import numpy as np
import sys             # For sys.argv, os.path.splitext and os.path.exist
import os


""" Extract HI Regions """
def extractHI(configOpts):
    # Extracts HI regions by calling the routine extractRegions (clumpfind?).
    # Skipped if region files already seem to exist.
    # Runs clumpfind on the HI regions? The regions are currently square.
    # Instead we call extractFromCoords.py with ds9RegsFile as input
    print
    print 'Extracting HI regions ... '
    hiImage = configOpts['hiImage']
    fuvImage = configOpts['fuvImage']
    ds9RegsFile = os.path.splitext(fuvImage)[0] + '.reg'
    # hiImage = "NGC_628_RO_MOM0_THINGS.FITS"
    hiMaskFile = configOpts['hiMaskFile'] # << but we want to extract based on UV sources instead
    #hiMaskFile = configOpts['fuvMaskFileRej']
    output = os.path.splitext(hiImage)[0]
    # Skip the extraction if the first file already exists
    #  !!hardcoded in extractRegions.py; hardcoded 001
    if os.path.exists('{0}_Regions/{0}_Reg001.fits'.format(output)):
        print "... extracted regions exist - skipping this step."
    else:
        if False in map(os.path.exists, [hiImage, hiMaskFile]):
            print 'Error: Either the image or the mask does not exist\n'
            sys.exit()
        print 'Warning: hardcoded HI region size of 20 pixels'
        pdr.extractFromCoords(hiImage, ds9RegsFile, output, imageDS9 = fuvImage, size=20, mosaic=False, mosaicBack=0)
    print

    #  Output: HI postage stamps and masks named as {output}_Regions/{output}_Msk###.fits
    return 0

""" FUV """
#FUV flux measurements
#First, try to read the flux file. Otherwise, calculate the fluxes. Output: fluxtable
def getflux(configOpts):

    print "Measuring FUV fluxes ..."
    uvheader, uvimage = pdr.open_image(configOpts['fuvImage'])
    data_uv = configOpts['data_uv']
    coords_uv = os.path.splitext(configOpts['fuvImage'])[0] + '_Peaks.dat' #as defined earlier

    try:
        dummy = file(data_uv, 'r')
        print "... fluxes read from file {0} instead.".format(data_uv)
        columns = np.dtype([('PDRID', 'int',), ('RA', '|S14'), ('DEC', '|S14'), ('aperture', 'float'), ('mean_at_r', 'float'), ('bgflux', 'float'), ('cumulflux', 'float'), ('netflux', 'float'), ('sflux', 'float')])
        fluxtable = pdr.read_array(data_uv, columns)
        print "  {0} records read".format(np.size(fluxtable.PDRID))
    except:
        #file does not exist, so we calculate. Force recalculation by removing the file.
        #Variables in: uvheader, uvimage, uvcoords
        #Variable out: fluxtable (includes uvcoords afterwards)
        #	fluxtable labels: PDRID, RA, DEC, aperture, mean_at_r, bgflux, cumulflux, netflux, sflux
        print "No previous flux measurements."
        uvcoords = pdr.read_coords(coords_uv)
        print "Flux will be measured at {0} positions on image {1} ...".format(len(uvcoords), configOpts['fuvImage'])
        print "Fluxes will be scaled by a factor {0}".format(configOpts['FUVscale'])
        #This scaling is useful since clumpfind doesn't seem to work properly with ~1e-15 values.
        fluxtable = pdr.fuv_flux(configOpts['fuvImage'], uvheader, uvcoords, configOpts, verbose = True)
        print "... done"

        #print fluxtable.RA, fluxtable.RA[0] #here fluxtable.RA[0] is simply a string

        #Write: UV fluxes
        fluxfile = fancyPrint(writeLog=True, logFile=data_uv, verbose=True)
        print "Flux info:"
        #print "#PDRID, aperture (arcsec), mean at r (units/arcsec^2), bg flux, cumul. flux, net flux"
        fluxfile.write("#PDRID, RA, DEC, aperture (arcsec), mean at r (units/pixel), bg flux, cumul. flux, net flux, sigma net flux\n")
        for n in range(len(fluxtable)):
            #fedora bug: changed PDRID:>3 to PDRID:3g
            #print "{0.PDRID:3g}, {1.aperture:5.2f}, {1.mean_at_r:7.5e}, {1.bgflux:7.5e}, {1.cumulflux:7.5e}, {1.netflux:7.5e}".format(uvcoords[n], fluxtable[n])
            fluxfile.write("{0.PDRID:3g}, {0.RA}, {0.DEC}, {1.aperture:7.5e}, {1.mean_at_r:7.5e}, {1.bgflux:7.5e}, {1.cumulflux:7.5e}, {1.netflux:7.5e}, {1.sflux:7.5e}\n".format(uvcoords[n], fluxtable[n]))
          #fluxfile.close()

    print
    #print "Exiting here."
    #exit(0)

    return uvheader, fluxtable

""" HI """
#HI patch measurements
#First, try to read the results file. Otherwise, detect patches with SExtractor (or Clumpfind).
#Output: hidata
def getHI(configOpts, fluxtable):
    print "Identifying HI patches ..."
    data_hi = configOpts['data_hi']
    hiImage = configOpts['hiImage']
    hi_baselist = "{0}_Regions/{0}_Reg".format(os.path.splitext(hiImage)[0])
    hi_logsbase = "{0}_Regions/{1}".format(os.path.splitext(hiImage)[0], configOpts['hi_logsbase'])
    try:
        dummy = file(data_hi, 'r')
        print "... HI patch data read from file {0}".format(data_hi)
        #Read: HI data (if already there); note that if the read fails, patches will be re-measured also
        columns = np.dtype([('PDRID', 'int',), ('RA', '|S14'), ('DEC', '|S14'), ('NHI', 'float'), ('sNHI', 'float')])
        hidata = pdr.read_array(data_hi, columns)
        print "  {0} records read".format(np.size(hidata))
        #Very simple and incomplete check for enough data:
        if (np.size(np.unique(hidata.PDRID)) != np.size(np.unique(fluxtable.PDRID))):
            print "Fatal error: number of UV and HI PDRIDs do not match. Try regenerating the HI files."
            exit(0)
    except:
        #File does not exist, so process HI postage stamp regions.
        #Must have regions with names matching the PDRIDs
        scale = configOpts['HIscale'] # 2.896e19 for NGC 628 (m?)Jy/beam to cm-2
        hi_bg = configOpts['hiBackground']
        print "Using map scaling of {0} cm-2".format(scale)
        hiresults = [[] for dummy in xrange(5)] #PDRID, ra, dec, NHI, sNHI
        #Warning: up to 999 files
        for i in range(np.size(fluxtable.PDRID)):
            fitsfile = hi_baselist + "{0:03d}.fits".format(int(fluxtable.PDRID[i]))
            print "Calling SExtractor with file {0}".format(fitsfile)

            catfile = pdr.call_SEx(fitsfile)
            if catfile==1:
                print "Fatal error: SExtractor call failed."
                exit(0)
            else:
                logfile = hi_logsbase + "{0:03}.txt".format(int(fluxtable.PDRID[i]))
                patches = pdr.read_secat(catfile, fitsfile, logfile, wcs=True)
                #Patches are not background-subtracted, that will happen now:

                for n in range(len(patches[0])):
                    if patches[2][n] > (hi_bg*2.):
                      #subtract the background and ignore patches that are fainter than hi_bg
                        hiresults[0].append(fluxtable.PDRID[i])
                        hiresults[1].append(patches[0][n])
                        hiresults[2].append(patches[1][n])
                        hiresults[3].append((patches[2][n]-hi_bg)*scale)
                        hiresults[4].append(hi_bg * scale * 0.5) #fix sN_HI to half the background value
                    else:
                        print "Patch too faint: ", patches[2][n]

                if not np.size(np.where(hiresults[0] == fluxtable.PDRID[i])):
                    print "No suitable patches were found for PDRID ", fluxtable.PDRID[i], "; creating empty entry."
                    #With SExtractor, this is unlikely since it does not pick up the faintest patches
                    hiresults[0].append(fluxtable.PDRID[i])
                    hiresults[1].append(0)
                    hiresults[2].append(0)
                    hiresults[3].append(0)
                    hiresults[4].append(0)

        #End for

        #Finally, convert the results into a rec array and save the results
        labels = 'PDRID, RA, DEC, NHI, sNHI'
        hidata = pdr.records(hiresults, labels)
        #Write: HI patches
        hifile = fancyPrint(writeLog=True, logFile=data_hi, verbose=True)
        print "HI patches:"
        #print "#PDRID, RA, DEC, NHI, sNHI"
        hifile.write("#PDRID, RA, DEC, NHI, sNHI\n")
        for n in range(len(hidata)):
            #print type(hidata[n].NHI)
            #Apparently python 2.6 needs explicit conversion from numpy.float32, so workaround
            #fedora bug: changed PDRID:>3 to PDRID:3g
            #print "{0.PDRID:3g}, {0.RA:14s}, {0.DEC:14s}, {1:7.5g}, {2:7.5g}".format(hidata[n], float(hidata[n].NHI), float(hidata[n].sNHI))
            hifile.write("{0.PDRID:3g}, {0.RA:14s}, {0.DEC:14s}, {1:7.5g}, {2:7.5g}\n".format(hidata[n], float(hidata[n].NHI), float(hidata[n].sNHI)))
          #hifile.close()
      #End of try/except

        print

        return hidata

""" RHO_HI & RGAL & G0 """
#Calculate rho_HI by combining UV and HI data. Also calculate Rgal in the same loop.
#Variables in: fluxtable, hidata, galpars, uvheader (for WCS)
#Variable out: rhoHI and srho (same length as HI data which could have multiple entries per UV entry); Rgal[kpc] (UV length), G0, contr (source contrast)
def getRhoRG(configOpts, fluxtable, hidata, uvheader):

    Rgal = []
    rhoHI = []
    srho = []
    n = 0
    c_RA =  configOpts['c_RA']
    c_DEC = configOpts['c_DEC']
    pa = configOpts['pa']
    incl = configOpts['incl']
    dist = configOpts['distance'] *1e6 #in Mpc but we want pc here
    #We estimate the error in rho to be of the order of half a pixel size
    pix_size = uvheader.getXPixelSizeDeg() * 3600 #Assume square pixel
    srho_fixed = dist / 3600. * np.pi / 180. * 0.5 * pix_size #srho in parsec

    print "Calculating Rho_HI, Rgal, G0"
    for i in range(np.size(fluxtable.PDRID)):
        Rgal.append(pdr.separation(fluxtable.RA[i], fluxtable.DEC[i], c_RA, c_DEC, pa, incl, dist*1e-3, uvheader)) #dist passed in kpc so Rgal will be in kpc
        while (n < np.size(hidata.PDRID)) and (hidata.PDRID[n] == fluxtable.PDRID[i]):
            rhoHI.append(pdr.separation(fluxtable.RA[i], fluxtable.DEC[i], hidata.RA[n], hidata.DEC[n], pa, incl, dist, uvheader)) #dist in pc
            srho.append(srho_fixed)
            n += 1

    #Variables in: fluxtable.netflux, fluxtable.mean_at_r, dist, pix_size, ext, and rho_HI
    #variable out: G0, contr
    #dist already extracted above; pix_size calculated above
    ext = configOpts['ext']
    G0 = []
    contr = [] #source contrast
    for i in range(np.size(hidata.PDRID)):
        uv_idx = np.where(fluxtable.PDRID == hidata.PDRID[i]) #one match expected only
        G0.append(pdr.Gnaught(fluxtable.netflux[uv_idx], dist, rhoHI[i], ext))
        contr.append(pdr.contrast(fluxtable.netflux[uv_idx], fluxtable.mean_at_r[uv_idx], rhoHI[i], dist))

    print
    #Now we have rhoHI, srho and G0 for each hidata.PDRID
    #G0 and contr end up as a list of arrays, that is unwieldy, so fix here:
    G0 = map(float, G0)
    contr = map(float, contr)

    return rhoHI, srho, Rgal, G0, contr

""" DUST TO GAS RATIO """
#Calculate dust and write to file (needs to be separated) - once per PDRID, that is once per UV source
#Variables in: fluxtable, hidata, Rgal; c_RA, c_DEC, dustmodel, pa, incl, dist
#Variables out: dd0, sdd0 (size of hidata, but at the RA, DEC of the fuvdata)
def dusttogas(configOpts, fluxtable, hidata, Rgal):

    print "Determining dust-to-gas ratios ..."
    c_RA =  configOpts['c_RA']
    c_DEC = configOpts['c_DEC']
    pa = configOpts['pa']
    incl = configOpts['incl']
    dist = configOpts['distance'] * 1e6 # [pc]

    dusttype, dustmodel = pdr.read_dust(configOpts)
      #dust parameters extracted from configOpts depending on the dust model (const, slope, map)

    dd0 = []; sdd0 = []
    for i in range(np.size(hidata.PDRID)):
        uv_idx = np.where(fluxtable.PDRID == hidata.PDRID[i])[0][0]
        #one match expected only -- but uv_idx is an array. So take [0][0], otherwise there will be issues using it
        #print uv_idx, fluxtable.RA[0], fluxtable.RA[uv_idx]#, fluxtable.RA[uv_idx[0]]
        #hack: fluxtable.RA[uv_idx] yields [' hh mm ss '] instead of plain string, which causes problems here. Hence [0]
        #dust routine expects D_gal(kpc) but ignored anyway?
        param_array = [[fluxtable.RA[uv_idx], fluxtable.DEC[uv_idx]], [c_RA, c_DEC], dustmodel, [pa, incl, dist*1e-3, Rgal[uv_idx]], uvheader]

        #print param_array
        d, s = pdr.determine_dust(dusttype, param_array) #! should separate dust model from calculating R
        if (s == 0.):
            print "Warning: using hardcoded values for sdd0." #Can add options of const or relative + value
            s = 0.1

        dd0.append(d)
        sdd0.append(s) #all errors are absolute

    print "Warning: using hardcoded 8.69 solar value and assuming 12+log(O/H) input for dust map"
    print

    #Now we have dd0, sdd0 for each hidata.PDRID

    return dd0, sdd0

""" RE-ORGANIZE AND LOG """
#Aggregate everything before calculating ntot
#Write and print the results by collating uv, hi data, and combined data, and create the 'data' table
def collate(configOpts, fluxtable, hidata, dd0, sdd0, rhoHI, srho, G0, contr):

    print "Re-organizing & tabulating ..."
    print "PDRID, RA, DEC, pRA, pDEC, NHI, sNHI, dd0, sdd0, rhoHI, srho, net flux, sflux, mean_at_r, ntot, sntot"
    data = [[] for dummy in range(20)]
    data_full = configOpts['data_full']
    ntot = np.zeros(np.size(hidata.PDRID))
    sntot = np.zeros(np.size(hidata.PDRID))
      #Add ntot and sntot for later
    with open(data_full, 'w') as fullfile: #no close file needed this way
        for i in range(np.size(hidata.PDRID)):
      #There can be no duplicate PDRIDs in fluxtable or this will yield incorrect results
          #for legibility:
          uv_idx = np.where(fluxtable.PDRID == hidata.PDRID[i])[0][0] #expect only one value - problem if there are more HI PDRIDs than UV PDRIDs
          #hack: 'g' does not take [number] so have to convert again with [0]. Need to 'fix' fluxtable. Same problem with G0
          #print "{0:>3}, {1:12s}, {2:12s}, {3:8g}, {4:6g}, {5:.3f}, {6:f}, {7:.3f}, {8:.3f}, {9:g}, {10:.3g}, {11:g}".format(hidata.PDRID[i], fluxtable.RA[uv_idx], fluxtable.DEC[uv_idx], hidata.NHI[i], hidata.sNHI[i], dd0[i], sdd0[i], rhoHI[i], srho[i], fluxtable.netflux[uv_idx], fluxtable.sflux[uv_idx], fluxtable.mean_at_r[uv_idx])

          #Patch RAs, DECs not written here
          #fullfile.write("{0:>3}, {1:12s}, {2:12s}, {3:8g}, {4:6g}, {5:.3f}, {6:f}, {7:.3f}, {8:.3f}, {9:g}, {10:.3g}, {11:g}\n".format(hidata.PDRID[i], fluxtable.RA[uv_idx], fluxtable.DEC[uv_idx], hidata.NHI[i], hidata.sNHI[i], dd0[i], sdd0[i], rhoHI[i], srho[i], fluxtable.netflux[uv_idx], fluxtable.sflux[uv_idx], fluxtable.mean_at_r[uv_idx]))
          #Fedora bug: need to type cast the variables
          fullfile.write("{0:>3}, {1:12s}, {2:12s}, {3:8g}, {4:6g}, {5:.3f}, {6:f}, {7:.3f}, {8:.3f}, {9:g}, {10:.3g}, {11:g}\n".format(
                int(hidata.PDRID[i]),
                fluxtable.RA[uv_idx], fluxtable.DEC[uv_idx],
                hidata.NHI[i], hidata.sNHI[i],
                dd0[i], sdd0[i], rhoHI[i], srho[i],
                fluxtable.netflux[uv_idx], fluxtable.sflux[uv_idx],
                fluxtable.mean_at_r[uv_idx]))

          for m, n in enumerate((hidata.PDRID[i], fluxtable.RA[uv_idx], fluxtable.DEC[uv_idx], hidata.RA[i], hidata.DEC[i], Rgal[uv_idx], hidata.NHI[i], hidata.sNHI[i], dd0[i], sdd0[i], rhoHI[i], srho[i], fluxtable.netflux[uv_idx], fluxtable.sflux[uv_idx], fluxtable.mean_at_r[uv_idx], fluxtable.aperture[uv_idx], G0[i], contr[i], ntot[i], sntot[i])):
              #print m, n
              data[m].append(n)

    labels = 'PDRID, RA, DEC, pRA, pDEC, Rgal, NHI, sNHI, dd0, sdd0, rhoHI, srho, flux, sflux, mean_at_r, aperture, G0, contr, ntot, sntot'
    data = pdr.records(data, labels) #transform into rec.array
    #print data #instead of the print earlier, for testing
    print "Wrote to file {0}".format(data_full)
    print "... done"
    print

    return data

""" CALCULATE NTOT """
#Final step: Calculate n, sigma_n
def calculatentot(configOpts, data):
    print "Calculating n_tot ..."

    #ext, dist already extracted from configOpts in rho_HI section
    R25 = configOpts['R25']
    ext = configOpts['ext']
    dist = configOpts['distance'] * 1e6 #[pc]

    for i in range(np.size(data.PDRID)):
        data.ntot[i] = pdr.PDRmodel(data.flux[i], ext, dist, data.dd0[i], data.rhoHI[i], data.NHI[i])
        data.sntot[i] = pdr.errors(data.ntot[i], data.flux[i], data.sflux[i], data.rhoHI[i], data.srho[i], data.dd0[i], data.sdd0[i], data.NHI[i], data.sNHI[i], ext, "errordump.log")

    return data

""" FILTER RESULTS"""
# These filters need to be carefully justified
def filter(data):
    pdb.set_trace()

    #Filter values here - a placeholder. Need to filter in justifiable ways before writing results.
    #print data.rhoHI
    print
    print "Results pre-filtering: ", np.size(data.PDRID)

    data = data.compress(data['rhoHI'] < 500) #plausibility cut-off
    #data = data.compress(data['rhoHI'] < 1000)
    print "Censoring rho_HI > 500 pc (plausilibity cut-off), remaining ", np.size(data.PDRID)

    data = data.compress(data['contr'] > 1)
    print "Censoring contrast < 1, remaining ", np.size(data.PDRID)

    data = data.compress(data['dd0'] > 0.001) #filter the fake values where no dd0 data available
    print "Filtering values without dd0 values, remaining ", np.size(data.PDRID)
    print
    return data

""" WRITE RESULTS """
def dump(configOpts, data):
    #For now, we will output n and coordinates
    #'Loose variables' coming in from earlier in the pipeline: G0 ###and Rgal
    resultcsv = fancyPrint(writeLog=True, logfile=configOpts['data_results'], verbose=True)
    R25 = configOpts['R25']
    print "Results: RA, DEC, Rgal, ntot, sntot, sntot/ntot (%)"
    for i in range(np.size(data.ntot)):
        #print "{0:12s}, {1:12s}, {2:7.2f}, {3:5.1e}, {4:5.1e}, {5:3.0f}".format(data.RA[i], data.DEC[i], data.Rgal[i], data.ntot[i], data.sntot[i], data.sntot[i] / data.ntot[i] *100)

        #Here we also print a csv for my plotting routines
        #The PDRID is last for compatibility reasons.
        #Again problem with 'g' and 'str' type in data:
        #print Rgal[i], Rgal[i]/galpars.R25, data.flux[i], data.rhoHI[i], data.NHI[i], data.G0[i], data.dd0[i], ntot[i], np.log10(ntot[i]), data.contr[i]
        #  resultcsv.write("{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g},{7:g},{8:g},{9:g}\n".format(data.Rgal[i], data.Rgal[i]/galpars.R25, data.flux[i], data.rhoHI[i], data.NHI[i], data.G0[i], data.dd0[i], ntot[i], np.log10(ntot[i]), data.contr[i], data.PDRID[i]))

        #Write the full results instead of the limited results above
        #'PDRID, RA, DEC, Rgal, NHI, sNHI, dd0, sdd0, rhoHI, srho, flux, sflux, mean_at_r, G0, contr', ntot, sntot
        #Fedora bug
        #resultcsv.write("{0.PDRID:>3},{0.RA:s},{0.DEC:s},{0.pRA:s}, {0.pDEC:s}, {0.Rgal:g},{1:g},{0.NHI:g},{0.sNHI:g},{0.dd0:g},{0.sdd0:g},{0.rhoHI:g},{0.srho:g},{0.flux:g},{0.sflux:g},{0.mean_at_r:g},{0.aperture:g},{0.G0:g},{0.contr:g},{0.ntot:g},{0.sntot:g}\n".format(data[i], data.Rgal[i]/R25))
        resultcsv.write("{0.PDRID:3g},{0.RA:s},{0.DEC:s},{0.pRA:s}, {0.pDEC:s}, {0.Rgal:g},{1:g},{0.NHI:g},{0.sNHI:g},{0.dd0:g},{0.sdd0:g},{0.rhoHI:g},{0.srho:g},{0.flux:g},{0.sflux:g},{0.mean_at_r:g},{0.aperture:g},{0.G0:g},{0.contr:g},{0.ntot:g},{0.sntot:g}\n".format(data[i], data.Rgal[i]/R25))


    print "Median error: ", np.median(data.sntot / data.ntot) * 100

    return 0

""" PRESENTATION """
def present(data):

    print "Writing DS9 region file with n"
    pdr.ntot_ds9(configOpts['densities_ds9'], data)
    print

    print "Creating a quick plot ..."

    R25 = configOpts['R25']
    pdr.plot_all(data.Rgal/R25, (data.dd0, data.flux, data.NHI, data.rhoHI, data.G0, data.ntot), r"$\rm{R/R_{25}}$", ("d/d0", "FUV Flux", "NHI (cm-2)", "rhoHI (pc)", "G0", "ntot"), (True, True, False, False, True, True))
    #Or just a R vs n plot
    #pdr.plot_n(data.Rgal/R25, np.log10(ntot))
    return 0


###############################################################################
###                                 MAIN                                    ###
###############################################################################

def main(configFile, createConfig=False, verbose=True, overwrite=False):
    """
    This main routine sequentially calls all the subroutines needed for the PDR
    method to work.
    """

    # configOptionsInstance is a instance of the class ConfigOptions. For convenience,
    # configOpts allows a quick access to the configOptionsInstance.options dict.
    configOptionsInstance = pdr.ConfigOptions(configFile, createConfig=createConfig,
                                              verbose=verbose)

    if createConfig is True:
        if configOptionsInstance.creationStatus is True:
            return
        else:
            raise IOError('Impossible to create %s' % configFile)
    else:
        configOpts = configOptionsInstance.options

    # Creates the logger

    if configOpts['log'] is True:
        logger = pdr.fancyPrint(logFile=configOpts['logFile'], writeLog=True,
                                verbose=verbose)
    else:
        logger = pdr.fancyPrint(logFile=None, writeLog=False, verbose=verbose)

    logger.write('Pipeline run with options:\n', doPrint=False)
    configOptionsInstance.logOpts(logger)

    # Adds the verbose and overwrite options to the configOpts dictionary
    configOpts['verbose'] = verbose
    configOpts['overwrite'] = overwrite

    # Calls the getRegions routine, which produces the region masks for the FUV and HI images
    pdr.getRegions(configOpts, logger)

    return

    #Todo: call reject regions on _Peaks.dat if desired

    extractHI(configOpts)

    uvheader, fluxtable = getflux(configOpts)

    hidata = getHI(configOpts, fluxtable)

    rhoHI, srho, Rgal, G0, contr = getRhoRG(configOpts, fluxtable, hidata, uvheader)

    dd0, sdd0 = dusttogas(configOpts, fluxtable, hidata, Rgal)

    data = collate(configOpts, fluxtable, hidata, dd0, sdd0, rhoHI, srho, G0, contr)

    data = calculatentot(configOpts, data)

    data = filter(data)

    dump(configOpts, data)

    present(data)  # This includes filtering values (a placeholder)

    return


if __name__ == '__main__':
    """
    Reads command line arguments
    """

    print

    from optparse import OptionParser

    usage = 'usage: python %prog configFile [runs pipeline]' + \
            '\n       python %prog -c configFile [generates a sample config file]'

    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--config', dest='sampleFile',
                      help='generates a default configuration file.',
                      default=False, action='store_true')
    parser.add_option('-q', '--quiet',
                      action='store_false', dest='verbose', default=True,
                      help='don\'t print status messages to stdout')
    parser.add_option('-o', '--overwrite',
                      action='store_true', dest='overwrite', default=False,
                      help='overwrites files. If not present, the pipeline tries ' +
                           'to skip the steps already completed.')

    (options, args) = parser.parse_args()

    if len(args) == 0:
        if options.sampleFile is False:
            parser.error('Incorrect number of arguments')
        else:
            main('PMAP.dat', createConfig=True, verbose=options.verbose)
    else:
        configFile = args[0]
        if options.sampleFile is True:
            main(configFile, createConfig=True, verbose=options.verbose)
        else:
            if os.path.exists(configFile):
                main(configFile, createConfig=False, verbose=options.verbose,
                     overwrite=options.overwrite)
            else:
                parser.error('File not found.')
