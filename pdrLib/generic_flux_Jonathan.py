def generic_flux(filename, x, y, maxrad_pix, pix_size = 1.0, contrast = 0.5, logfile = False, verbose = False):
  from calcFlux import calcFlux

  #Calculate flux in concentric rings (for each coordinate) (qphot can do this)
  fluxlist = []
  avglist = []
  mean_at_r = 0; aperture = 0; bgflux = 0; cumulflux = 0; netflux = 0; sflux = 0
  for n in range(maxrad_pix): #iterate over rings of width 1 pixel
    r = n + 1
    fluxlist.append(calcFlux(filename, x, y, r))
    
    if n > 0:
      ring_avg = (fluxlist[n][1] - fluxlist[n-1][1]) / (fluxlist[n][0] - fluxlist[n-1][0]) #incremental flux divided by ring surface area
      avglist.append(ring_avg) #current avglist index will be n-1 since it starts at n=1
      if (logfile):
        logdump(logfile, "{0:3d} {1:7.2e} {2:8.3e} {3:7.5e}\n".format(r, fluxlist[n][0], fluxlist[n][1], ring_avg))

    if (n > 1) and (not mean_at_r):
      #First check if the mean value has dropped to at least 50% (internals.contrast) of the initial value (arbitrary)
      #  note that if the actual max value is much higher, 50% is not indicative of any kind of source contrast
      #Then check if a local minimum was encountered.
      #This could problems with very poor source contrast, but those may need to be filtered out earlier on
      #For example: need to drop to at least half the initial mean value first
      if (avglist[n-1] < (avglist[0] * contrast)):
        if (avglist[n-1] > avglist[n-2]): #list is 1 shorter than n
          mean_at_r = avglist[n-2] #per pixel
          aperture = (r-1) * pix_size #in arcsec and assuming square pixels
          bgflux = mean_at_r * fluxlist[n-2][0] #flux/pixel*area[in pixels]
          cumulflux = fluxlist[n-2][1]
          netflux = cumulflux - bgflux
          mean_at_r = mean_at_r / (pix_size*pix_size) #!need to report per arcsec^2
          sflux = 0.3*netflux #assign an uncertainty on the measurement here !!hardcoded for now!!
          if (logfile):
            logdump(logfile, "#First local minimum registered\n")
          #print aperture, r-1

  if (not mean_at_r) or (verbose): #we exit leaving all the values at zero
    print "Warning: minimum not found!"
    if (logfile): 
      logdump(logfile, "#Warning: minimum not found!\n")

  return mean_at_r, aperture, bgflux, cumulflux, netflux, sflux
