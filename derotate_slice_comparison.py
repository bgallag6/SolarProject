# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 21:41:33 2017

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as patches

Z = np.load('C:/Users/Brendan/Desktop/derotate_mpi_4seg.npy')

# define subregion coordinates   
x1 = -550
x2 = 550
y1 = -500
y2 = 500

coords_type = 'arc'

# create a list of all the files. This is USER-DEFINED
flist = glob.glob('F:/Users/Brendan/Desktop/SolarProject/FITS2/aia*.fits')
nf = len(flist)

# Select the image that is the "middle" of our selection.
# We do this because the solar derotation algorithm operates centered on the 
mid_file = np.int(np.floor(nf / 2))

mc_list = []  # create an empty list

count = 0  # counter to see where program is at 

# Use defined coordinates, extract the submaps from each AIA image, and store
# them in the empty list. This takes many minutes to complete.
print " "
print "Reading files and extracting submaps. This takes a while..."
print " "
   
if coords_type == 'arc':
    for filename in flist:
        mc_list.append(Map(filename).submap([x1,x2]*u.arcsec, [y1,y2]*u.arcsec))
        if (count%10) == 0:
           print("file %i out of %i" % (count,nf))  # just to see where program is at
        count += 1
    
	
# Create a new Map Cube object to hold our de-rotated data
new_mapcube = Map(mc_list, cube=True)
print " "
print "Creating derotated cube..."
print "Please wait..."

# Perform the derotation of the submaps. This take a while too.
#dr = mapcube_solar_derotate(new_mapcube)
dr = mapcube_solar_derotate(new_mapcube, layer_index=mid_file)

mid_subarr = dr[mid_file].data		# extract data from middle file of derotated datacube
rem_i = mid_subarr.shape[0]  # calculate remainder to get integer dimensions
rem_j = mid_subarr.shape[1]  # calculate remainder to get integer dimensions
subarr_idim = mid_subarr.shape[0]
subarr_jdim = mid_subarr.shape[1]
 
mc_list = np.zeros((1,2,3))  # free up memory
new_mapcube = np.zeros((1,2,3))  # free up memory

t = dr[0].date  # extract the date / time from the first image
base_time = (t.hour * 60.) + (t.minute * 60.) + t.second  # convert date / time to seconds

# initialize arrays to hold exposure time, pixel data, and time values
I = np.empty((nf))  # exposure time
DATA = np.empty((nf, subarr_idim, subarr_jdim), dtype=np.int16)  # save as int16, since that is what original is
TIME = np.empty((nf))  # might as well just have all as float

# loop through datacube and extract pixel data and time values
for p in range(0,nf):
    Ex = dr[p].exposure_time
    I[p] = Ex.value
    L = dr[p].data
    L_trim = L[0:mid_subarr.shape[0], 0:mid_subarr.shape[1]]
    DATA[p][:][:] = L_trim  # normalize by exposure time
    T = dr[p].date
    curr_time=(T.hour * 3600.)+(T.minute * 60.)+T.second	
    TIME[p] = curr_time - base_time  # calculate running time of image
    


for i in range(DATA.shape[0]):
    fig = plt.figure(figsize=(20,10))
    ax1 = plt.subplot2grid((1,11),(0, 0), colspan=5, rowspan=1)
    ax1.set_title('1-Segment Derotation', fontsize=20, y=1.01)
    rect = patches.Rectangle((1150,0), 50, 1678, color='white', fill=False)  
    ax1.add_patch(rect)
    rect2 = patches.Rectangle((350,0), 50, 1678, color='white', fill=False)  
    ax1.add_patch(rect2)
    plt.imshow(np.flipud(DATA[i]), cmap='sdoaia171', vmin=25, vmax=1000)
    plt.xticks(visible=False)
    ax2 = plt.subplot2grid((1,11),(0, 6), colspan=5, rowspan=1)
    ax2.set_title('4-Segment Derotation', fontsize=20, y=1.01)
    rect3 = patches.Rectangle((1150,0), 50, 1678, color='white', fill=False)  
    ax2.add_patch(rect3)
    rect4 = patches.Rectangle((350,0), 50, 1678, color='white', fill=False)  
    ax2.add_patch(rect4)
    plt.imshow(np.flipud(Z[i]), cmap='sdoaia171', vmin=25, vmax=1000)
    plt.savefig('C:/Users/Brendan/Desktop/derotate_mpi_4/hour_%i.jpeg' % i)
    
    