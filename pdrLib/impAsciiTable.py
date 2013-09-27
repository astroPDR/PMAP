
from Error import raiseError

try:
    import asciitable as at
except:
    try:
        import astropy.io.ascii as at
    except:
        raiseError('No asciitable module found.', None)
