# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 14:24:12 2018

@author: Brendan
"""
import os
import yaml


stream = open('specFit_config.yaml', 'r')
cfg = yaml.load(stream)

directory = cfg['fits_dir']
date = cfg['date']
wavelength = cfg['wavelength']


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
