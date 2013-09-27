# SExtractor functions

import subprocess  # to call SExtractor
# Ferguson's code to read a SExtractor catalog
from sextutils import se_catalog
from . import pf, pw
import numpy as np

# Run SExtractor, extract the results
# SExtractor requires files default.conv default.nnw pdrSExout.param and HIclumps.sex
# pdrLib requires only that pdrSEx.param include ALPHAPEAK_J2000 and
# DELTAPEAK_J2000

# Call SExtractor - this procedure is simply a specialized wrapper to SExtractor
# Need SExtractor in path
# It needs to be called for every postage stamp HI region


def call_SEx(fitsfile):
    # The HIclumps.sex config file comes with this software. It does not need to be customized further.
    # Output catalog: HIclumps.cat
    err = subprocess.Popen(["sex", fitsfile, '-c',  'HIclumps.sex'], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE).wait()

    if (err):
        # print "Error calling SExtractor."
        return 1
    # If no clumps were found, the .cat file will only contain commented lines
    # This file only needs to be hard-coded here and in HIclumps.sex
    return "HIclumps.cat"

# Call Fergusons's library to read a SExtractor catalog,
# produced previously to select HI patches in an individual candidate PDRs


def read_secat(catfile, fitsfile):
    patchtable = [[] for dummy in xrange(3)]  # ra, dec, NHI
    # Must catch empty catalog file:
    # Call will fail if the catalog file is empty or does not have the
    # expected format
    try:
        catalog = se_catalog(catfile, readfile=True, preserve_case=False)
    except:
        # print "Warning: no clumps found. Creating empty entry."
        return None

    records = len(catalog.alphapeak_j2000)
    hduImage = pf.open(fitsfile)
    wcs = pw.WCS(hduImage[0].header)
    for n in range(records):
        a, d = catalog.alphapeak_j2000[n], catalog.deltapeak_j2000[n]
        x, y = wcs.wcs_world2pix(np.array([[a, d]]), 0)[0]
        x = int(np.round(x))
        y = int(np.round(y))
        maxHI = hduImage[0].data[y, x]
        # no offset? Verify!!
        # print "a,d: ", a, d, "x,y: ", x, y, "value: ", maxHI
        patchtable[0].append(a)
        patchtable[1].append(d)
        patchtable[2].append(maxHI)
        # patchtable[3].append(0.1)    #error now calculated after calling this routine
        # print "Warning: hardcoding sNHI at 0.1"
        # logdump(logfile, "{} {} {}\n".format(a, d, maxHI))

    return patchtable  # calling function administrates PDRID
