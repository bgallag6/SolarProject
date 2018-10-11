# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 19:58:41 2017

@author: Brendan
"""

import os
import sys

# set variables from command line
directory = sys.argv[1]
date = sys.argv[2]
wavelength = sys.argv[3]
#directory = 'S:'
#date = '20101213'
#wavelength = 1700

"""
# make FITS directory
outdir = '%s/FITS/%s/%i/' % (directory, date, wavelength)

if not os.path.exists(os.path.dirname(outdir)):
    try:
        print("Making output directory structure...")
        os.makedirs(os.path.dirname(outdir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST: raise
"""        
    
# make FITS directory
outdir = '%s/FITS/%s/' % (directory, date)

if wavelength == 'UV':
    if not os.path.exists(os.path.dirname(outdir+'1600/')):
        try:
            print("Making output directory structure...")
            os.makedirs(os.path.dirname(outdir+'1600/'))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST: raise
    if not os.path.exists(os.path.dirname(outdir+'1700/')):
        try:
            print("Making output directory structure...")
            os.makedirs(os.path.dirname(outdir+'1700/'))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST: raise
            
elif wavelength in ['171', '193', '211', '304', '1600', '1700', 'magnetogram', 'continuum']:
    if not os.path.exists(os.path.dirname(outdir+wavelength+'/')):
        try:
            print("Making output directory structure...")
            os.makedirs(os.path.dirname(outdir+wavelength+'/'))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST: raise