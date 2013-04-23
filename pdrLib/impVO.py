try:
    from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
    voModule = True
except:
    try:
        from vo.tree import VOTableFile, Resource, Table, Field
        voModule = True
    except:
        voModule = False
