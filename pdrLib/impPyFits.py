
from Error import raiseError

try:
    import pyfits as pf
except:
    try:
        import astropy.io.fits as pf
    except:
        raiseError('No pyfits module found.', None)
