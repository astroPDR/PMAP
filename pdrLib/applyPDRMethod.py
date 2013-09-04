#!/usr/bin/env python
# encoding: utf-8
"""
applyPDRMethod.py

Created by Jose Sanchez-Gallego on 30 Aug 2013.
Copyright (c) 2013. All rights reserved.

This file compiles the functions relative to the PDR method calculations.
Originally written by Jonathan Heiner. Compiled and updated on 08/29/2013
by José Sánchez-Gallego.

"""

from . import pf, pw
import numpy as np
from calcSeparation import calcSeparation
from . import Gnaught, contrast, read_dust, determine_dust, PDRmodel, errors
from . import table
from Error import raiseWarning
import os


def getRhoRG(configOpts, fluxFUVTable, hiData, logger):
    """
    Calculate rho_HI by combining UV and HI data. Also calculate Rgal in the same loop.
    Variables in: fluxFUVTable, hiData, galpars
    Variable out: rhoHI and srho (same length as HI data which could have
    multiple entries per UV entry); Rgal[kpc] (UV length), G0, contr (source
    contrast)
    """

    dist = configOpts['distance']
    header = pf.getheader(configOpts['fuvImage'])
    wcsFUV = pw.WCS(header)
    try:
        pix_size = np.array([np.abs(wcsFUV.wcs.cd[0, 0]), np.abs(wcsFUV.wcs.cd[1, 1])]) * 3600
    except:
        pix_size = np.abs(wcsFUV.wcs.cdelt[0:2]) * 3600

    dataRho = table.Table(None, names=('PDRID', 'HIID', 'rhoHI', 'sRhoHI', 'Rgal', 'G0', 'Contrast'),
                          dtype=('i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8'))

    c_RA = configOpts['c_RA']
    c_DEC = configOpts['c_DEC']
    pa = configOpts['pa']
    incl = configOpts['incl']
    dist = configOpts['distance'] * 1e6  # in Mpc but we want pc here
    ext = configOpts['ext']
    # We estimate the error in rho to be of the order of half a pixel size
    srho_fixed = dist / 3600. * np.pi / 180. * 0.5 * pix_size.mean()  # srho in parsec

    logger.write('Calculating Rho_HI, Rgal, G0 ...', newLine=True)

    for ii in range(len(fluxFUVTable)):

        pdrID = fluxFUVTable['PDRID'][ii]

        Rgal = calcSeparation(
            fluxFUVTable[ii]['RA'], fluxFUVTable[ii]['Dec'], c_RA, c_DEC,
            pa, incl, dist * 1e-3,  # dist passed in kpc so Rgal will be in kpc
            wcsFUV)

        hiPatchesPDR = hiData[hiData['PDRID'] == pdrID]
        for hiRow in hiPatchesPDR:
            hiID = hiRow['HIID']
            rhoHI = calcSeparation(
                fluxFUVTable[ii]['RA'], fluxFUVTable[ii]['Dec'],
                hiRow['RA'], hiRow['Dec'],
                pa, incl, dist, wcsFUV)  # dist in pc

            G0 = Gnaught(fluxFUVTable[ii]['netflux'], dist, rhoHI, ext)
            contr = contrast(fluxFUVTable[ii]['netflux'], fluxFUVTable[ii]['mean_at_r'],
                             rhoHI, dist)

            dataRho.add_row([pdrID, hiID, rhoHI, srho_fixed, Rgal, G0, contr])

    return dataRho


def dustToGas(configOpts, fluxFUVTable, hiData, dataRho, logger):
    """
    DUST TO GAS RATIO
    Calculate dust and write to file (needs to be separated) - once per PDRID,
    that is once per UV source.
    Variables in: fluxFUVTable, hiData, Rgal; c_RA, c_DEC, dustmodel, pa, incl, dist
    Variables out: dd0, sdd0 (size of hiData, but at the RA, DEC of the fuvdata)
    """

    logger.write('Determining dust-to-gas ratios ...', newLine=True)

    dataDust = dataRho.copy()

    header = pf.getheader(configOpts['fuvImage'])
    wcsFUV = pw.WCS(header)

    c_RA = configOpts['c_RA']
    c_DEC = configOpts['c_DEC']
    pa = configOpts['pa']
    incl = configOpts['incl']
    dist = configOpts['distance'] * 1e6  # [pc]

    dusttype, dustmodel = read_dust(configOpts)
    # dust parameters extracted from configOpts depending on the dust model (const, slope, map)

    dd0 = []
    sdd0 = []
    warningIssued = False
    for rowHI in hiData:
        pdrID = rowHI['PDRID']
        # print uv_idx, fluxFUVTable.RA[0], fluxFUVTable.RA[uv_idx]#, fluxFUVTable.RA[uv_idx[0]]

        # hack: fluxFUVTable.RA[uv_idx] yields [' hh mm ss '] instead
        # of plain string, which causes problems here. Hence [0]
        # dust routine expects D_gal(kpc) but ignored anyway?

        pdrTable = fluxFUVTable[fluxFUVTable['PDRID'] == pdrID]
        Rgal = dataRho[dataRho['PDRID'] == pdrID]['Rgal'][0]

        param_array = [[pdrTable['RA'], pdrTable['Dec']], [c_RA, c_DEC],
                       dustmodel, [pa, incl, dist * 1e-3,
                                   Rgal, wcsFUV]]

        # print param_array
        #! should separate dust model from calculating R
        d, s = determine_dust(dusttype, param_array)
        if s == 0. and warningIssued is False:
            # Can add options of const or relative + value
            raiseWarning('Warning: using hardcoded values for sdd0.', logger, newLine=False)
            s = 0.1
            warningIssued = True

        dd0.append(d)
        sdd0.append(s)  # all errors are absolute

    raiseWarning(
        'Warning: using hardcoded 8.69 solar value and assuming 12+log(O/H) input for dust map',
        logger, newLine=False)

    # Now we have dd0, sdd0 for each hiData.PDRID
    tableDD0 = table.Table([dd0, sdd0], names=('dd0', 'sdd0'))
    dataDust.add_columns(tableDD0.columns.values())

    return dataDust


def collate(configOpts, fluxFUVTable, hiData, dataDust, logger):
    """
    Compiles all the data from different tables and creates a single one. Two
    additional fields are created (empty) for NTot and SNTot
    """

    logger.write('Compiling all data ... ', newLine=True)

    data = dataDust.copy()

    # Adds coordinate columns for both FUV and HI
    coordTable = table.Table(np.zeros((len(data), 4)),
                             names=('FUV_RA', 'FUV_Dec', 'HI_RA', 'HI_Dec'),
                             dtype=4*['S20'])
    data.add_columns(coordTable.columns.values(), indexes=[2, 2, 2, 2])

    # Concatenates the data from hiData
    data = table.join(data, hiData['PDRID', 'HIID', 'NHI', 'sNHI'], keys=('PDRID', 'HIID'))

    # Adds columns for the FUV data and NTot
    fuvCols = table.Table(np.zeros((len(data), 6)),
                          names=('aperture', 'mean_at_r', 'netflux', 'sflux', 'NTot', 'sNTot'),
                          dtype=6*['f8'])
    data.add_columns(fuvCols.columns.values())

    # Now, lets complete the empty fields
    for row in data:
        fuvRow = fluxFUVTable[fluxFUVTable['PDRID'] == row['PDRID']]
        hiRow = hiData[(hiData['PDRID'] == row['PDRID']) & (hiData['HIID'] == row['HIID'])]

        row['FUV_RA'] = fuvRow['RA']
        row['FUV_Dec'] = fuvRow['Dec']
        row['HI_RA'] = hiRow['RA']
        row['HI_Dec'] = hiRow['Dec']

        for field in ['aperture', 'mean_at_r', 'netflux', 'sflux']:
            row[field] = fuvRow[field]

    # Renames netflux to flux
    data['netflux'].name = 'flux'

    return data


def calculateNTot(configOpts, data, logger):
    """
    Final step: Calculate n, sigma_n
    """

    logger.write('Calculating N(total) ...')

    # ext, dist already extracted from configOpts in rho_HI section
    ext = configOpts['ext']
    dist = configOpts['distance'] * 1e6  # [pc]

    for ii in range(len(data)):
        data['NTot'][ii] = PDRmodel(
            data['flux'][ii], ext, dist, data['dd0'][ii], data['rhoHI'][ii], data['NHI'][ii])
        data['sNTot'][ii] = errors(
            data['NTot'][ii], data['flux'][ii], data['sflux'][ii], data['rhoHI'][ii], data['sRhoHI'][ii],
            data['dd0'][ii], data['sdd0'][ii], data['NHI'][ii], data['sNHI'][ii], ext, logfile=False)

    return data


def filterData(data, logger):
    """
    FILTER RESULTS
    These filters need to be carefully justified
    """

    logger.write('Filtering results ... ', newLine=True)
    # Filter values here - a placeholder. Need to filter in justifiable ways before writing results.

    logger.write('Results pre-filtering: %d' % len(data))

    # Removes HI patches beyond a certain point
    filteredData = data[data['rhoHI'] < 500]  # plausibility cut-off
    # filteredData = data[data['RhoHI'] < 1000]

    logger.write('Censoring rho_HI > 500 pc (plausilibity cut-off), remaining %d' % len(filteredData))

    filteredData = filteredData[filteredData['Contrast'] > 1]
    logger.write('Censoring contrast < 1, remaining %d' % len(filteredData))

    # Filter the fake values where no dd0 data available
    filteredData = filteredData[filteredData['dd0'] > 0.001]
    logger.write('Filtering measurements without dd0 values, remaining %d' % len(filteredData))

    return filteredData


def applyPDRMethod(configOpts, fluxFUVTable, hiData, logger):

    dataRho = getRhoRG(configOpts, fluxFUVTable, hiData, logger)

    # Check rhoHIs, seem too large (>1kpc). The problem may be in the use of wcs conversions in the
    # calcSeparation routine. It should not be necessary.

    dataDust = dustToGas(configOpts, fluxFUVTable, hiData, dataRho, logger)

    # Creates the final table with all the compiled data and two empty columns where NTot and SNTot will go.
    data = collate(configOpts, fluxFUVTable, hiData, dataDust, logger)

    # We calculate N(total)
    data = calculateNTot(configOpts, data, logger)

    filteredData = filterData(data, logger)

    # Saves the data
    logger.write('Saving data ...', newLine=True)
    data.write(configOpts['data_results_unfiltered'], format='ascii.fixed_width', delimiter=' ')
    if configOpts['save_unfiltered_vot'] is True:
        fileVOT = os.path.splitext(configOpts['data_results_unfiltered'])[0] + '.vot'
        if os.path.exists(fileVOT):
            os.remove(fileVOT)
        data.write(fileVOT, format='votable')

    filteredData.write(configOpts['data_results'], format='ascii.fixed_width', delimiter=' ')
    if configOpts['save_results_vot'] is True:
        fileVOT = os.path.splitext(configOpts['data_results'])[0] + '.vot'
        if os.path.exists(fileVOT):
            os.remove(fileVOT)
        filteredData.write(fileVOT, format='votable')

    return data
