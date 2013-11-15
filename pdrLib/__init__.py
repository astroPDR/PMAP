
from Error import raiseError

try:
    import pyfits as pf
except:
    try:
        import astropy.io.fits as pf
    except:
        raiseError('No pyfits module found.', None)

try:
    import pywcs as pw
except:
    try:
        import astropy.wcs as pw
    except:
        raiseError('No pywcs module found.', None)

try:
    import asciitable as at
except:
    try:
        import astropy.io.ascii as at
    except:
        raiseError('No asciitable module found.', None)

try:
    import astropy.table as table
except:
    raiseError('No asciitable module found.', None)

from pdrLib import *
from getRegions import *
from createCatalog import *
from parseOptions import *
from fancyPrint import *
from fluxes import *
from getHIFlux import *
from applyPDRMethod import applyPDRMethod
from calcSeparation import calcSeparation
from calcScale import calcScale
from matchFOV import matchFOV, maskImage
