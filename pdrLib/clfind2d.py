#!/usr/bin/env python
# encoding: utf-8
"""
clfind2d.py

Created by José Ramón Sánchez-Gallego on 2010-06-07.

This script is a translation of the clfind2d.pro code by J. P. Williams
www.ifa.hawaii.edu/~jpw/clfind/.

Needs: pyfits, numpy and the routine search2d provided in a separated file
Usage: clfind2d.py fileIn levels [-l|--nolog] [-o|--out]
        fileIn: The original FITS file
        levels: The levels of the contours [l1,l2,l3,...]
        -l|--nolog: Disable the logging system
        -o|--out: The root of the output files. If not present, the [fileIn] root
                 will be used
        -n|-npixmin: Minimum number of pixels for a region to be accepted
        -i|--interval: Interval of levels [start,end,step]. Disables the levels argument

"""


import numpy as np
import pyfits as pf
import sys
import os
from time import strftime, gmtime, clock
from scipy.ndimage import generate_binary_structure, label


def printLog(text, log, verbose=True):
    if verbose is True:
        sys.stdout.write(text)
        sys.stdout.flush()
    if log is not None:
        unit = open(log, 'a')
        print >>unit, text,
        unit.close()
    return


# This routine identifies all the clumps between two consecutive contour levels.
def defReg(data, levMin, levMax, diagonal=False):

# Pixels with values between the two levels
    levelMask = (data >= levMin) & (data < levMax)
    nPixInLevel = np.sum(levelMask)

# Is no pixels between [levMin,levMax) are found, returns
    if nPixInLevel == 0:
        return nPixInLevel, None, 0

    if diagonal:
        struct = generate_binary_structure(2, 2)
    else:
        struct = generate_binary_structure(2, 1)

    reg, nReg = label(levelMask, structure=struct)

    return nPixInLevel, reg, nReg


# Decides if a set of pixels is a new clump or belong to a previous clump.
# In that case, merges the two clumps or mask the pixels to the closest
# clump
def defClump(data, levMin, reg, nReg, mask, clumpPeak, nCl):

    nNew = 0

    if nReg == 0: return 0

    nPix2, pix2, nReg2 = defReg(data, levMin, 99999, diagonal=True)

    for nR in range(1, nReg + 1):

        pix1 = np.where(reg == nR)
        nPix1 = pix1[0].size

        peak = np.argmax(data[pix1])
        iPeak = pix1[0][peak]
        jPeak = pix1[1][peak]

        maskPix2 = pix2 == pix2[iPeak, jPeak]
        aPix2 = mask[maskPix2]

        if np.amax(aPix2) == 0:
            # Found new clump!
            nNew += 1
            mask[maskPix2] = nCl + nNew
            peak = np.argmax(data[pix1])
            iPeak = pix1[0][peak]
            jPeak = pix1[1][peak]
            clumpPeak[nCl+nNew, 0] = iPeak
            clumpPeak[nCl+nNew, 1] = jPeak
        else:
            nC = np.unique(np.sort(aPix2, axis=None))
            if np.amin(nC) == 0:
                if nC.size == 2:
                    # Extending clump
                    mask[maskPix2] = nC[1]
                else:
                    clumpMerge = nC[1:]
                    # Merging clumps
                    iC = clumpPeak[clumpMerge, 0]
                    jC = clumpPeak[clumpMerge, 1]
                    for nR1 in range(0, nPix1):
                        i = pix1[0][nR1]
                        j = pix1[1][nR1]
                        d = (i-iC)**2 + (j-jC)**2
                        m = np.argmin(d)
                        mask[i, j] = clumpMerge[m]

    return nNew


# Clumps with fewer pixels than nPixMin are removed
# The labels of the clumps are resorted accordingly
def testBad(data, nPixMin, nCl, clumpPeak, mask):

    dMax = data[clumpPeak[1:nCl+1, 0], clumpPeak[1:nCl+1, 1]]
    newOrder = np.argsort(dMax)[::-1]

    nClNew = 0
    nBad = 0
    mask0 = mask.copy()

    for n0 in newOrder:
        iClp = mask0 == n0+1
        count = np.sum(iClp)

        if count <= nPixMin:
            nBad += 1
            mask[iClp] = 0
        else:
            nClNew += 1
            mask[iClp] = nClNew

    return nBad


# Main routine
def clfind2d(file, mask, levels, log=True, nPixMin=20, verbose=True):

# The final name of the mask
    fileOut = mask

# If log == True, sets the name of the log file
    if log is True:
        logFile = fileOut + '.log'
        if os.path.exists(logFile): os.remove(logFile)
    else:
        logFile = None

    tStart = clock()  # The clockwatch is on
    currentTime = strftime('%a, %d %b %Y %H:%M:%S GMT', gmtime())

    printLog('\n----------------------------------------------------------------\n', logFile, verbose=verbose)
    printLog('CLFIND2d: %s\n' % currentTime, logFile, verbose=verbose)
    printLog('----------------------------------------------------------------\n', logFile, verbose=verbose)
    printLog('Filename = %s\n' % file, logFile, verbose=verbose)
    printLog('Mask file = %s\n' % fileOut, logFile, verbose=verbose)
    if log is True:
        printLog('Log file = %s\n' % logFile, logFile, verbose=verbose)

    printLog('Minimum number of pixels = %d\n' % nPixMin, logFile, verbose=verbose)
    levStr = np.round(levels, 4)
    levStr = ','.join(map(str, levStr))
    printLog('Contour levels at = %s\n' % levStr, logFile, verbose=verbose)

    printLog('----------------------------------------------------------------\n', logFile, verbose=verbose)

# Loads the FITS image
    data = pf.getdata(file)
    header = pf.getheader(file)

    iSize, jSize = data.shape

# Bad pixels are set to -999.9
    data[np.isnan(data)] = -999.9

# Adds a top level to the contours list
    levels = np.append(levels, 99999)
    nLev = levels.size

    max_clumps = 1e6

# The array which will contain the mask
    mask = np.zeros([iSize, jSize], np.int)
    clumpPeak = np.zeros([max_clumps, 2], np.int)

# For each contour level, identifies the clumps and sorts them
    nCl = 0
    for nWork in range(nLev-2, -1, -1):

        levMin = levels[nWork]
        levMax = levels[nWork + 1]

        nPix, reg, nReg = defReg(data, levMin, levMax, diagonal=True)

        printLog('Contour level %f: %i pixels, %i regions, ' % (levels[nWork], nPix, nReg), logFile, verbose=verbose)

        nNew = defClump(data, levels[nWork], reg, nReg, mask, clumpPeak, nCl)

        printLog('%i new clumps\n' % (nNew), logFile, verbose=verbose)
        nCl += nNew

# Tests for bad pixels
    nBad = testBad(data, nPixMin, nCl, clumpPeak, mask)
    nCl -= nBad

    printLog('%d clumps found (%d rejected)\n' % (nCl, nBad), logFile, verbose=verbose)
    printLog('================================================================\n', logFile, verbose=verbose)

# Saves the mask
    if os.path.exists(fileOut): os.remove(fileOut)
    hdu = pf.PrimaryHDU(mask)
    hdu.header = header
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(fileOut)

    printLog('Writing output file: %s\n' % fileOut, logFile, verbose=verbose)

    tEnd = clock()
    interval = (tEnd - tStart)/60.0
    printLog('%.1f minutes elapsed\n' % interval, logFile, verbose=verbose)
    # printLog('\n', logFile, verbose=verbose)

    return


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: %prog [options] fileIn levels"
    parser = OptionParser(usage=usage)
    parser.add_option('-o', '--out', dest='out',
                      help='root of the output files')
    parser.add_option('-l', '--nolog', dest='log', action='store_false',
                      help='doesn\'t write the log file', default=True)
    parser.add_option('-n', '--npixmin', dest='npix',
                      help='minimum number of pixels for a region to be accepted', type=int, default=20)
    parser.add_option('-i', '--interval', dest='interval',
                      help='interval of levels [start,end,step]')

    (options, args) = parser.parse_args()

    if len(args) != 2 and not options.interval:
        parser.error('Incorrect number of arguments')

    fileIn = args[0]

    if not options.interval:
        try:
            levels = np.array(map(np.float32, args[1].split(',')))
        except:
            parser.error('Levels not understood')
    else:
        try:
            interv = map(np.float, options.interval.split(','))
            levels = np.arange(interv[0], interv[1] + interv[2] * 0.5, interv[2])
        except:
            parser.error('Interval not understood')

    if options.out:
        rootOut = options.out
    else:
        rootOut = os.path.splitext(os.path.basename(fileIn))[0]

    clfind2d(fileIn, rootOut, levels, log=options.log, nPixMin=options.npix)
