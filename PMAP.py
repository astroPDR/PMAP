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
import os


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
    if 'overwrite' not in configOpts:
        configOpts['overwrite'] = overwrite
    else:
        if configOpts['overwrite'] is not True:
            configOpts['overwrite'] = overwrite
        else:
            pass

    # Calls the getRegions routine, which produces the region masks for the FUV and HI images
    pdr.getRegions(configOpts, logger)

    #Todo: call reject regions on _Peaks.dat if desired

    fluxFUVTable = pdr.getUVFlux(configOpts, logger)
    hiData = pdr.getHIFlux(configOpts, fluxFUVTable, logger)

    # Now we apply the PDR mehod to the FUV and HI measurements
    pdrData = pdr.applyPDRMethod(configOpts, fluxFUVTable, hiData, logger)

    return pdrData


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
