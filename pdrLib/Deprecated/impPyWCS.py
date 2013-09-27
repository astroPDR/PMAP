
from Error import raiseError

try:
    import pywcs as pw
except:
    try:
        import astropy.wcs as pw
    except:
        raiseError('No pywcs module found.', None)
