# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 11:46:12 2018

@author: Brendan
"""

"""
#######################
#######################
# Generates derotated datacube of subregion.
#######################
#######################
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import glob
import sunpy
from sunpy.map import Map
from sunpy.physics.solar_rotation import mapcube_solar_derotate
import astropy.units as u
from astropy.time import Time
import yaml


stream = open('specFit_config.yaml', 'r')
cfg = yaml.load(stream)

directory = cfg['fits_dir']
date = cfg['date']
wavelength = cfg['wavelength']
sub_reg_coords = cfg['sub_reg_coords']
coords_type = cfg['coords_type']
bin_frac = cfg['bin_frac']
mmap_derotate = cfg['mmap_derotate']
save_temp = cfg['save_temp']

    
# rebin region to desired fraction 
def rebin(a, *args):
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    return eval(''.join(evList))
    
# define subregion coordinates   
x1 = sub_reg_coords[0]  # j1
x2 = sub_reg_coords[1]  # j2
y1 = sub_reg_coords[2]  # i1
y2 = sub_reg_coords[3]  # i2

# create a list of all the files. This is USER-DEFINED
flist = sorted(glob.glob('%s/FITS/%s/%i/aia*.fits' % (directory,date,wavelength)))
#flist = sorted(glob.glob('S:/FITS/%s/%i/aia*.fits' % (date,wavelength)))
nf = len(flist)

# Select the image that is the "middle" of our selection.
# We do this because the solar derotation algorithm operates centered on the 
mid_file = np.int(np.floor(nf / 2))

mc_list = []  # create an empty list

count = 0  # counter to see where program is at 

c3 = SkyCoord((x1)*u.arcsec, y1*u.arcsec, frame=frames.Helioprojective)  
c4 = SkyCoord((x2)*u.arcsec, y2*u.arcsec, frame=frames.Helioprojective) 

# Use defined coordinates, extract the submaps from each AIA image, and store
# them in the empty list. This takes many minutes to complete.
print(" ", flush=True)
print("Reading files and extracting submaps. This takes a while...", flush=True)
print(" ", flush=True)
   
if coords_type == 'pix':
    for filename in flist:
        mc_list.append(Map(filename).submap([x1,x2]*u.pixel, [y1,y2]*u.pixel))
        if (count%10) == 0:
            print("file %i out of %i" % (count,nf), flush=True)  # just to see where program is at
        count += 1

if coords_type == 'arc':
    for filename in flist:
        mc_list.append(Map(filename).submap(c3,c4))
        if (count%10) == 0:
           print("file %i out of %i" % (count,nf), flush=True)  # just to see where program is at
        count += 1
    
	
# Create a new Map Cube object to hold our de-rotated data
new_mapcube = Map(mc_list, cube=True)
print(" ", flush=True)
print("Creating derotated cube...", flush=True)
print("Please wait...", flush=True)

# Perform the derotation of the submaps. This take a while too.
dr = mapcube_solar_derotate(new_mapcube, layer_index=mid_file)  # derotate around middle file

print("done derotating", flush=True)


mid_subarr = dr[mid_file].data		# extract data from middle file of derotated datacube
rem_i = mid_subarr.shape[0] % bin_frac  # calculate remainder to get integer dimensions
rem_j = mid_subarr.shape[1] % bin_frac  # calculate remainder to get integer dimensions
subarr_idim = (mid_subarr.shape[0] - rem_i) // bin_frac  # get ydim
subarr_jdim = (mid_subarr.shape[1] - rem_j) // bin_frac  # get xdim
 
del mc_list # free up memory
del new_mapcube # free up memory

# initialize arrays to hold exposure time, pixel data, and time values
I = np.empty((nf))  # exposure time
DATA = np.empty((nf, subarr_idim, subarr_jdim), dtype=np.int16)  # save as int16, since that is what original is
TIME = np.empty((nf))  # time stamps
AVG = np.zeros((subarr_idim, subarr_jdim))

# loop through datacube and extract pixel data and time values
# this is probably another loop that I could take out and extract directly
for p in range(nf):
    Ex = dr[p].exposure_time
    I[p] = Ex.value
    L = dr[p].data
    L_trim = L[0:(mid_subarr.shape[0] - rem_i), 0:(mid_subarr.shape[1] - rem_j)]
    #small_L = rebin(L_trim, L_trim.shape[0]//bin_frac, L_trim.shape[1]//bin_frac)
    #DATA[p][:][:] = small_L
    DATA[p][:][:] = L_trim
    AVG += (L_trim / Ex.value)  # create normalized average visual image
    T = dr[p].date
    TIME[p] = Time(T).jd  # extract julian day time from each image
    
TIME -= TIME[0]  # set all values equal to time since first entry
TIME = np.around(TIME*86400)  # get the time value in seconds, and round to nearest whole number

# save the time-array, and exposure-time datacubes as numpy files
np.save('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength), TIME)
np.save('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength), I)

# calculate the average-intensity image of the timeseries 
AVG /= nf

# determine the middle file of the timeseries
mid_num = (DATA.shape[0]//2)
mid = DATA[mid_num]

# store average and middle images in array
visual = np.zeros((2,AVG.shape[0],AVG.shape[1]))
visual[0] = AVG
visual[1] = mid

# save averaged visual image
np.save('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength), visual)

# generate images of each visual region, to see if as expected
#fig = plt.figure(figsize=(20,20))
#plt.imshow(visual[0])
#fig = plt.figure(figsize=(20,20))
#plt.imshow(visual[1])


if mmap_derotate == "y":
    orig_shape = np.array([DATA.shape[0], DATA.shape[1], DATA.shape[2]])
    
    # create memory-mapped array with similar datatype and shape to original array
    mmap_arr = np.memmap('%s/DATA/Temp/%s/%i/derotated_mmap.npy' % (directory, date, wavelength), dtype='%s' % DATA.dtype, mode='w+', shape=tuple(orig_shape))
    
    # write data to memory-mapped array
    mmap_arr[:] = DATA[:]
    
    # save memory-mapped array dimensions to use when loading
    np.save('%s/DATA/Temp/%s/%i/derotated_mmap_shape.npy' % (directory, date, wavelength), orig_shape)

    # save original array if specified
    if save_temp == "y":
        np.save('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength), DATA)
    
    # flush memory changes to disk, then remove memory-mapped object and original array
    del mmap_arr
    del DATA
    
else:
    np.save('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength), DATA)