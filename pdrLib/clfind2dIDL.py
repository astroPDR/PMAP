#!/usr/bin/env python
# encoding: utf-8
"""
clfind2dIDL.py

Created by José Ramón Sánchez-Gallego on 2010-06-11.

This script is a frontend for the clfind2d.pro IDL script. It creates a 
batch file, launches IDL and renames the resulting files.
    
"""

# The main routine
def clfind2dIDL(file, levels, nPixMin=20, log=True, verbose=False):
    
    import numpy as np
    import sys, os
    import subprocess
    import time
    import tempfile
    
    batchFile = tempfile.NamedTemporaryFile(delete=False)
    
# Creates the batch file with the instructions which will be executed in IDL    
    # unitBatch = open(batchFile, 'w')
    print >>batchFile, '.rnew %s/clfind2d' % os.path.dirname(os.path.realpath(__file__))
    levStr = ','.join(map(str, levels.tolist()))
    print >>batchFile, 'clfind2d,file=\'%s\',levels=[%s],npixmin=%i,/log' % (os.path.splitext(file)[0], levStr, nPixMin)
    print >>batchFile, 'exit'
    batchFile.close()

# Runs the batch file in IDL
    stdoutPipe = sys.stdout
    if verbose == False:
        stdoutPipe = open('/dev/null', 'w')
        sys.stdout.write('Running IDL ... ')
        sys.stdout.flush()
    
    try:
        cmd = 'idl {0}'.format(batchFile.name)
        run = subprocess.Popen(cmd, stdout=stdoutPipe, stderr=subprocess.STDOUT, shell=True)
        returnCode = run.wait()    
    except:
        print 'problem found'
        return [False,]
        
    if verbose == False: print 'done'
    
# Renames the files

    fileMskIDL = os.path.splitext(file)[0] + '.fits.clf'
    fileLogIDL = 'clfind2d.log'
        
    return [True, fileMskIDL, fileLogIDL]
    

