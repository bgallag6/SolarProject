# -*- coding: utf-8 -*-
"""
Created on Wed Dec 07 22:19:52 2016

@author: Brendan
"""

# want to include this in module - in datacube creation - so can get full subregion desired

import scipy.signal
import matplotlib
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from matplotlib.widgets import  RectangleSelector
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
from scipy.interpolate import interp1d
from scipy import signal
import numpy as np
import scipy.misc
import astropy.units as u
import h5py
from sunpy.physics.transforms.solar_rotation import calculate_solar_rotate_shift


# create a list of all the files. This is USER-DEFINED

flist = glob.glob('F:/Users/Brendan/Desktop/SolarProject/FITS/20120923/171b/aia*.fits')

# Create an empty list
mc_list = []

# Use coordinates we just defined, extract the submaps from each AIA image, and store
# then in the empty list. This takes many minutes to complete.
print " "
print "Reading files and extracting submaps. This takes a while..."
print " "
for filename in flist:
   mc_list.append(Map(filename))  # ** SEE NOTE BELOW ***
	
# NOTE: the u.pixel means we're selecting data based on pixel coordinates. Alternate coordinates
#		would be degrees of solar latitude/longitude positions
	
# Create a new Map Cube object to hold our de-rotated data
new_mapcube = Map(mc_list, cube=True)


shifts = calculate_solar_rotate_shift(new_mapcube, layer_index=1)

print shifts