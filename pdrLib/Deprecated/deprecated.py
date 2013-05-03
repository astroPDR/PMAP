# Default options dictionary.
defOpts = dict(distance=None,
               log=True,
               logFile='test.log',
               verbose=False,
               fuvImage=None,
               hiImage=None,
               outputdir='./',
               fuvMaskFile='Default',
               fuvInteractive=True,
               fuvSigma=None,
               fuvhLevel=None,
               fuvnLevels=None,
               checkFUV=True,
               minDistance=500,
               scale=None,
               joinRegions=True,
               fluxContrast=3.0,
               peakContrast=3.0,
               alwaysRemoveClose=True,
               fuvMaskFileRej='Default',
               hiMaskFile='Default',
               hiInteractive=True,
               hiSigma=None,
               hihLevel=None,
               hinLevels=None)


# Opens the configuration file and retrieved the values. Returns a dictionary
def readConfigFile(fileOpts):
    configDic = {}
    lines = open(fileOpts, 'r').read().splitlines()
    for line in lines:
        strippedLine = line.strip()
        if strippedLine == '' or strippedLine[0] == '#': continue
        opt = strippedLine.split()[0]
        value = ' '.join(strippedLine.split()[1:])
        if '#' in value:
            comment = value[value.find('#'):].strip()
            value = value[0:value.find('#')].strip()
        if value == '':
            value = None
        elif value.lower()  in ['y', 'true', 'yes']:
            value = True
        elif value.lower()  in ['n', 'no', 'false']:
            value = False
        elif value.lower() in ['none', '', '-']:
            value = None
        else:
            try:
                value = float(value)
            except:
                pass
        configDic[opt] = value
    return configDic


# Parses the command-line call and returns the FUV and HI image names and the config file
# if any. In that case, it reads the config file. Otherwise, loads the default options.
def parseCommandLine():

    global configOpts

    import optparse

    parser = optparse.OptionParser()

    usage = 'usage: %prog [options] [FUV_Image] [HI_Image]'
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-c', '--config', dest='configFile',
        help='configuration file', default=None)
    parser.add_option('-y', '--default', dest='default',
        action='store_true', help='does not prompt. Uses default answers',
        default=False)

    (options, args) = parser.parse_args()

# Checks if the config file exists and loads it. Otherwise loads the default dir
    if options.configFile != None:
        if not os.path.exists(options.configFile):
            fancyPrint('Config file not found. Using defaults.', newLine=True)
            configFile = None
            configOpts = defOpts.copy()
        else:
            configFile = options.configFile
            configOpts = readConfigFile(configFile)
            fancyPrint('Configuration file %s loaded' % configFile,
                newLine=True)
    else:
        fancyPrint('Loading defaults options', newLine=True)
        configOpts = defOpts.copy()

    if os.path.exists(configOpts['logFile']): os.remove(configOpts['logFile'])

# The command-line images take precedence over the ones in the config file (in case
# they exist)
    if len(args) >= 1 and args[0] != 'None':
        configOpts['fuvImage'] = args[0]
    if len(args) >= 2 and args[1] != 'None':
        configOpts['hiImage'] = args[1]

# Checks for potential errors in the input
    if configOpts['fuvImage'] == None or configOpts['fuvImage'] == None:
        fancyPrint('FUV image or HI image not defined\n', newLine=True)
        sys.exit()

    if not os.path.exists(configOpts['fuvImage']) or not os.path.exists(configOpts['fuvImage']):
        fancyPrint('Some of the files don\'t exist\n', newLine=True)
        sys.exit()

    configOpts['default'] = options.default

    return configFile


# As parseCommandLine, but in this case it parses the call when getRegions is called as
# a subroutine
def parseInputData(inputData):

    global configOpts

    # Copies the input dictionary
    extraData = inputData.copy()

    # If the configuration file is defined in the inputData dictionary, it
    # loads it. Otherwise, uses the default options
    if 'configFile' in inputData:
        configFile = inputData['configFile']
        del extraData['configFile']
        try:
            configOpts = readConfigFile(configFile)
        except:
            fancyPrint('Config file not found. Using defaults.', newLine=True)
            configFile = None
            configOpts = defOpts.copy()
    else:
        configFile = None
        configOpts = defOpts.copy()

    # Sets the default option to false
    configOpts['default'] = False

    # Loads the remaining options in inputData
    for key,value in extraData.items():
        configOpts[key] = value

    # Checks for potential errors in the input
    if configOpts['fuvImage'] == None or configOpts['fuvImage'] == None:
        fancyPrint('FUV image or HI image not defined\n', newLine=True)
        sys.exit()

    if not os.path.exists(configOpts['fuvImage']) or not os.path.exists(configOpts['fuvImage']):
        fancyPrint('Some of the files don\'t exist\n', newLine=True)
        sys.exit()

    return configFile


