# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 05:19:28 2017

@author: Brendan
"""

from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
import numpy as np
import astropy.units as u
from mpi4py import MPI



def datacube(directory, date, wavelength, sub_reg_coords, coords_type, bin_frac):
    
    # rebin region to desired fraction 
    def rebin(a, *args):
        shape = a.shape
        lenShape = len(shape)
        factor = np.asarray(shape)/np.asarray(args)
        evList = ['a.reshape('] + \
                 ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
                 [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
        #print ''.join(evList)
        return eval(''.join(evList))
        
    # define subregion coordinates   
    x1 = sub_reg_coords[0]  # j1
    x2 = sub_reg_coords[1]  # j2
    y1 = sub_reg_coords[2]  # i1
    y2 = sub_reg_coords[3]  # i2
    
    # create a list of all the files. This is USER-DEFINED
    flist = glob.glob('%s/FITS/%s/%i/aia*.fits' % (directory,date,wavelength))
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
   
    if coords_type == 'pix':
        for filename in flist:
            mc_list.append(Map(filename).submap([x1,x2]*u.pixel, [y1,y2]*u.pixel))
            if (count%10) == 0:
                print("file %i out of %i" % (count,nf))  # just to see where program is at
            count += 1
    
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
    
    print "done derotating"
    
    
    mid_subarr = dr[mid_file].data		# extract data from middle file of derotated datacube
    rem_i = mid_subarr.shape[0] % bin_frac  # calculate remainder to get integer dimensions
    rem_j = mid_subarr.shape[1] % bin_frac  # calculate remainder to get integer dimensions
    subarr_idim = (mid_subarr.shape[0] - rem_i) / bin_frac  # get ydim
    subarr_jdim = (mid_subarr.shape[1] - rem_j) / bin_frac  # get xdim
     
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
        L_trim = L[0:(mid_subarr.shape[0] - rem_i), 0:(mid_subarr.shape[1] - rem_j)]
        small_L = rebin(L_trim, L_trim.shape[0]/bin_frac, L_trim.shape[1]/bin_frac)
        DATA[p][:][:] = small_L  # normalize by exposure time
        T = dr[p].date
        curr_time=(T.hour * 3600.)+(T.minute * 60.)+T.second	
        TIME[p] = curr_time - base_time  # calculate running time of image
        
    return DATA, I, TIME
 
   
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

import sys

directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])
xmin = int(sys.argv[4])
xmax = int(sys.argv[5])
ymin = int(sys.argv[6])
ymax = int(sys.argv[7])
#coords_type = sys.argv[5]
coords_type = 'arc'
#bin_frac = int(sys.argv[6])
bin_frac = 1 


y_range = [i for i in range(ymin,ymax)]

chunks = np.array_split(y_range, size)  # Split the data based on no. of processors

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == 0:
        sub_reg_coords = [xmin, xmax, chunks[0][0], chunks[0][len(chunks[i])-1]+10]
    elif rank == (size-1):
        sub_reg_coords = [xmin, xmax, chunks[size-1][0]-10, chunks[size-1][len(chunks[i])-1]]
    elif rank == i:
        sub_reg_coords = [xmin, xmax, chunks[i][0]-10, chunks[i][len(chunks[i])-1]+10]
print sub_reg_coords       
"""
data_cube, exp_arr, time_arr = datacube(directory, date, wavelength, sub_reg_coords, coords_type, bin_frac)
derotated_cube = comm.gather(data_cube, root=0) # Gather all the results

# Again, just have one node do the last bit
if rank == 0:
  stack_cube = np.vstack(derotated_cube)
  print stack_cube.shape			# Verify we have a summed version of the input cube
  print exp_arr.shape
  print time_arr.shape
"""
"""    
# save the pixel-value, time-array, and exposure-time datacubes as numpy files
#np.save('%s/DATA/Temp/%s/%i/%i_%ii_%i_%ij_data_rebin%i.npy' % (directory, date, wavelength, y1, y2, x1, x2, bin_frac), DATA)
#np.save('%s/DATA/Temp/%s/%i/%i_%ii_%i_%ij_time.npy' % (directory, date, wavelength, y1, y2, x1, x2), TIME)
#np.save('%s/DATA/Temp/%s/%i/%i_%ii_%i_%ij_exposure.npy' % (directory, date, wavelength, y1, y2, x1, x2), I)
np.save('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength), DATA)
np.save('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength), TIME)
np.save('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength), I)

# calculate the average-intensity image of the timeseries 
AVG = np.average(DATA,axis=0)

# determine the middle file of the timeseries
mid_num = (DATA.shape[0]/2)
mid = DATA[mid_num]

print "Middle file is number %i" % mid_num

print "(100,100): %i ~ %i" % (AVG[100][100], mid[100][100])  # check values are reasonably close

# check the two image sizes agree
print " Average Image Dimensions = %i, %i" % (AVG.shape[0], AVG.shape[1])
print " Middle Image Dimensions = %i, %i" % (mid.shape[0], mid.shape[1])

# store average and middle images in array
visual = np.zeros((2,AVG.shape[0],AVG.shape[1]))
visual[0] = AVG
visual[1] = mid

print visual.shape  # check array size agrees with expected

# save visual-image array
#np.save('%s/DATA/Output/%s/%i/%i_%ii_%i_%ij_visual.npy' % (directory, date, wavelength, y1, y2, x1, x2), visual)
np.save('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength), visual)

# generate images of each visual region, to see if as expected
#fig = plt.figure(figsize=(20,20))
#plt.imshow(visual[0])
#fig = plt.figure(figsize=(20,20))
#plt.imshow(visual[1])
"""

