# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 15:10:30 2017

@author: Brendan
"""


import os
import sys

# set variables from command line
directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])
#directory = 'S:'
#date = '20101213'
#wavelength = 1700

"""
# make FITS directory
outdir = '%s/FITS/%s/%i/' % (directory, date, wavelength)

if not os.path.exists(os.path.dirname(outdir)):
    try:
        print "Making output directory structure..."
        os.makedirs(os.path.dirname(outdir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST: raise
"""

# make temp directory
outdir = '%s/DATA/Temp/%s/%i/' % (directory, date, wavelength)

if not os.path.exists(os.path.dirname(outdir)):
    try:
        print("Making output directory structure...")
        os.makedirs(os.path.dirname(outdir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST: raise

# make output directory
outdir = '%s/DATA/Output/%s/%i/' % (directory, date, wavelength)

if not os.path.exists(os.path.dirname(outdir)):
    try:
        print("Making output directory structure...")
        os.makedirs(os.path.dirname(outdir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST: raise
            