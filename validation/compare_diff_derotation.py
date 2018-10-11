# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 08:06:45 2018

@author: Brendan
"""

import glob
import sunpy
from sunpy.map import Map
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
from sunpy.physics.solar_rotation import calculate_solar_rotate_shift
from sunpy.physics.differential_rotation import diffrot_map
import matplotlib.pyplot as plt
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate

from timeit import default_timer as timer
from mpi4py import MPI



"""
directory = 'S:'
date = '20130626'
wavelength = 1600
#x1 = -200
#x2 = 200
#y1 = -500
#y2 = -100
x1 = -500
x2 = 500
y1 = -500
y2 = 500

flist = sorted(glob.glob('%s/FITS/%s/%i_shift/aia*.fits' % (directory,date,wavelength)))

nf = len(flist)

# Select the middle image, to derotate around
mid_file = np.int(np.floor(nf / 2))

# make mapcube containing first and last maps & calculate derotation shifts
# if not centered at 0deg long, shift amount wont be enough -- maybe only calculate latitude, use other method for longitude trim
mc_shifts = []
mapI = Map(flist[0])
mapF = Map(flist[-1])

mc_shifts.append(mapI)
mc_shifts.append(mapF)
new_mapcube1 = Map(mc_shifts, cube=True)
shifts = calculate_solar_rotate_shift(new_mapcube1, layer_index=0)
diffrot_longitude = np.abs(np.floor((np.floor(shifts['x'][1].value)/2.)))  # since reflects first-last frame shift, divide by 2 for shifts around center frame
diff_long_pix = diffrot_longitude / (mapI.scale)[0].value  # calculate rotation amount in pixels
diffrot_lat_pix = int((np.abs((shifts['y'][1].value)) / (mapI.scale)[1].value) * 5.)  # calculate total latitude shift in pixels, x2 since underestimated?

c3 = SkyCoord((x1-diffrot_longitude)*u.arcsec, y1*u.arcsec, frame=frames.Helioprojective)  
c4 = SkyCoord((x2+diffrot_longitude)*u.arcsec, y2*u.arcsec, frame=frames.Helioprojective) 

# get middle frame subregion & time to anchor derotation
midmap = Map(flist[mid_file]).submap(c3,c4)
dt0 = midmap.date

# calculate pixels to trim based off of what actual derotation trims *for some reason this result is different than method above
diff_mapI = diffrot_map(Map(flist[0]).submap(c3, c4), time=dt0)
diff_mapF = diffrot_map(Map(flist[-1]).submap(c3, c4), time=dt0)

xminindI = np.argmin(np.fliplr(diff_mapI.data),axis=1)[diffrot_lat_pix:-diffrot_lat_pix]
xminindF = np.argmin(diff_mapF.data,axis=1)[diffrot_lat_pix:-diffrot_lat_pix]

xminI = midmap.data.shape[1] - np.min(xminindI)  
xminF = midmap.data.shape[1] - np.min(xminindF)

dCube = np.empty((nf, midmap.data.shape[0], midmap.data.shape[1]), dtype=np.int16)  # save as int16, since that is what original is

count = 0

# loop through datacube and extract pixel data and time values
for filename in flist:
    smap = Map(filename).submap(c3, c4)
    #dmap = diffrot_map(smap, dt=-dt*u.second)
    dmap = diffrot_map(smap, time=dt0)
    dCube[count] = dmap.data
    count += 1        
    

dCube_trim = dCube[:,diffrot_lat_pix:-diffrot_lat_pix,xminI:-xminF]



    
mc_list = []  # create an empty list

count = 0  # counter to see where program is at 

for filename in flist:
    mc_list.append(Map(filename).submap(c3,c4))
    count += 1
    	
# Create a new Map Cube object to hold our de-rotated data
new_mapcube = Map(mc_list, cube=True)

# Perform the derotation of the submaps. This take a while too.
dr = mapcube_solar_derotate(new_mapcube, layer_index=mid_file)

mid_subarr = dr[mid_file].data		# extract data from middle file of derotated datacube
rem_i = mid_subarr.shape[0] % 1  # calculate remainder to get integer dimensions
rem_j = mid_subarr.shape[1] % 1  # calculate remainder to get integer dimensions
subarr_idim = int((mid_subarr.shape[0] - rem_i) / 1)  # get ydim
subarr_jdim = int((mid_subarr.shape[1] - rem_j) / 1)  # get xdim
 
mc_list = np.zeros((1,2,3))  # free up memory
new_mapcube = np.zeros((1,2,3))  # free up memory

# initialize arrays to hold exposure time, pixel data, and time values
data = np.empty((nf, subarr_idim, subarr_jdim), dtype=np.int16)  # save as int16, since that is what original is

for p in range(nf):
    L = dr[p].data
    L_trim = L[0:(mid_subarr.shape[0] - rem_i), 0:(mid_subarr.shape[1] - rem_j)]
    data[p][:][:] = L_trim
"""

#"""
#dc3 = dCube_trim[-1] - dCube_trim[0]

directory = 'S:'
date = '20110909'
wavelength = 1700

cube_shape = np.load('%s/DATA/Temp/%s/%i/derotated_mmap_shape.npy' % (directory, date, wavelength))
cube = np.memmap('%s/DATA/Temp/%s/%i/derotated_mmap.npy' % (directory, date, wavelength), dtype='int16', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))

cube_shape_non = np.load('%s/DATA/Temp/%s/%i/derotated_mmap_shape_non.npy' % (directory, date, wavelength))
cube_non = np.memmap('%s/DATA/Temp/%s/%i/derotated_mmap_non.npy' % (directory, date, wavelength), dtype='int16', mode='r', shape=(cube_shape_non[0], cube_shape_non[1], cube_shape_non[2]))

#dc = np.load('S:/DATA/Temp/20110909/1700/derotated.npy')
#dc2 = np.load('S:/DATA/Temp/20140822/1600/derotated.npy')  

for i in range(12):
    fig = plt.figure(figsize=(18,10))
    (ax1, ax2) = fig.subplots(1, 2)
    
    """
    ax1.imshow(np.flipud(data[i,15:-5]), cmap='sdoaia1600', vmin=50, vmax=500)
    ax1.vlines(820,0,1615,'r')
    ax1.hlines(1275,0,1625,'r')
    ax2.imshow(np.flipud(dCube_trim[i]), cmap='sdoaia1600', vmin=50, vmax=500)
    ax2.vlines(820,0,1620,'r')
    ax2.hlines(1275,0,1668,'r')
    plt.savefig('C:/Users/Brendan/Desktop/derot/Bframe%i.jpeg' % i, bbox_inches='tight')
    """
    
    ax1.imshow(np.flipud(cube[i*50][:,30:-37]), cmap='sdoaia1700', vmin=500, vmax=3000)
    ax1.vlines(420,0,1200,'r')
    ax1.hlines(375,0,1150,'r')
    ax1.vlines(820,0,1200,'r')
    ax1.hlines(875,0,1150,'r')
    ax2.imshow(np.flipud(cube_non[i*50][2:-2]), cmap='sdoaia1700', vmin=500, vmax=3000)
    ax2.vlines(420,0,1215,'r')
    ax2.hlines(375,0,1150,'r')
    ax2.vlines(820,0,1215,'r')
    ax2.hlines(875,0,1150,'r')
    #plt.savefig('C:/Users/Brendan/Desktop/derot/Cframe%i.jpeg' % i, bbox_inches='tight')

#plt.imshow(np.flipud(dc3))
#"""
