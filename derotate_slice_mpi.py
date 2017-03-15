# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 00:06:41 2017

@author: Brendan
"""

from timeit import default_timer as timer

import matplotlib.pyplot as plt
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.solar_rotation import mapcube_solar_derotate
import numpy as np
from astropy import units as u
#import astropy.units as u
from sunpy.physics.solar_rotation import calculate_solar_rotate_shift   # fix these module imports in other scripts as well - get rid of 'transform'
from mpi4py import MPI



def datacube(directory, date, wavelength, sub_reg_coords, coords_type, bin_frac):
    """
    Generates derotated datacube of subregion.
    
    directory : 
        Path to folder containing .fits files - also where cube is saved.  (String)
        
    date :
        Date of timeseries [YYYY/MM/DD]  (String)
        
    wavelength :
        Wavelength of dataset.  (Int)
        
    sub_reg_coords : 
        Coordinates of subregion [x1,x2,y1,y2] (int/float)
    
    coords_type :
        Arcsecond or Pixel coordinates ['arc' or 'pix'] (String)
                   
    bin_frac : 
        The fraction by which the image should be rebinned by. (int)
      
    Example:
    ::
        ss.datacube(directory='F:/SDO/data/20130530/1600', date='20130530', wavelength=1600,
           sub_reg_coords=[2200,3000,2300,2600], coords_type='pix', bin_frac=2) 
    """
    
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
        
    x1a = sub_reg_coords[0] # j1
    x2a = sub_reg_coords[1]  # j2
    y1a = sub_reg_coords[2]  # i1
    y2a = sub_reg_coords[3]  # i2
    
    # Create an empty list
    mc_shifts = []
    #flist_shift = glob.glob('%s/FITS/%s/%i/aia*.fits' % (directory,date,wavelength))
    #flist_shift = glob.glob('F:/Users/Brendan/Desktop/SolarProject/derotate_test/aia*.fits') # works
    #flist_shift = glob.glob('F:/Users/Brendan/Desktop/SolarProject/FITS2/aia*.fits')
    flist = glob.glob('%s/FITS/%s/%i/aia*.fits' % (directory,date,wavelength))
    #flist = glob.glob('%s/aia*.fits' % (directory))
    #nf1 = len(flist_shift)
    nf = len(flist)
    print nf
    # Select the image that is the "middle" of our selection.
    # We do this because the solar derotation algorithm operates centered on the 
    mid_file = np.int(np.floor(nf / 2))
    
    # Select the image that is the "middle" of our selection.
    # We do this because the solar derotation algorithm operates centered on the 
    #mid_file1 = np.int(np.floor(nf1 / 2))
    
    if coords_type == 'pix':
        #for filename in flist_shift:
        for filename in flist:
            if filename == flist[0] or filename == flist[mid_file] or filename == flist[nf-1]:
                mc_shifts.append(Map(filename).submap([x1a,x2a]*u.pixel, [y1a,y2a]*u.pixel))

    
    if coords_type == 'arc':
        #for filename in flist_shift:
        for filename in flist:
            if filename == flist[0] or filename == flist[mid_file] or filename == flist[nf-1]:
                mc_shifts.append(Map(filename).submap([x1a,x2a]*u.arcsec, [y1a,y2a]*u.arcsec))

            
    #print nf1, mid_file1
    
    new_mapcube1 = Map(mc_shifts, cube=True)
    
    shifts = calculate_solar_rotate_shift(new_mapcube1, layer_index=1)
    M = shifts['x']
    L = shifts['y']
    
    # Create a new Map Cube object to hold our de-rotated data
    
    
    
    
    print sub_reg_coords   
    # define subregion coordinates   
    x1 = sub_reg_coords[0]-np.array(M[0])  # j1
    x2 = sub_reg_coords[1]-np.array(M[len(M)-1])  # j2
    y1 = sub_reg_coords[2]-np.array(L[0])-10  # i1
    y2 = sub_reg_coords[3]-np.array(L[len(M)-1])+10  # i2
    print x1,x2,y1,y2
    
    # create a list of all the files. This is USER-DEFINED
    #flist = glob.glob('%s/FITS/%s/%i/aia*.fits' % (directory,date,wavelength))
    #flist = glob.glob('%s/aia*.fits' % (directory))
    #nf = len(flist)
    #print nf
    # Select the image that is the "middle" of our selection.
    # We do this because the solar derotation algorithm operates centered on the 
    #mid_file = np.int(np.floor(nf / 2))
    
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
    print len(mc_list)
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
        
    print DATA.shape
        
    # generate images of each visual region, to see if as expected
    #fig = plt.figure(figsize=(20,20))
    #plt.imshow(AVG)
    return DATA, TIME, I

"""
def spec_fit(coords):
    x1 = coords[0]
    x2 = coords[1]
    y1 = coords[2]
    y2 = coords[3]
    print x1,x2,y1,y2        
"""
     
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
	
start = timer()

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

#directory = 'F:/Users/Brendan/Desktop/SolarProject/FITS2'
#date = '20130626'
#wavelength = 171

import sys

directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])

x_min = int(sys.argv[4])
x_max = int(sys.argv[5])

y_min = int(sys.argv[6])
y_max = int(sys.argv[7])



#x_coords = [-500,500]
x_coords = [x_min, x_max]
y_coords = [i for i in range(y_min,y_max)]


chunks = np.array_split(y_coords, size)  # Split the data based on no. of processors

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]
ymin = subcube[0]
ymax = subcube[len(subcube)-1]
#print ymin, ymax
coords = [x_coords[0],x_coords[1],ymin,ymax]
coords_type = 'arc'
bin_frac = 1
#testr = datacube(directory, date, wavelength, coords, coords_type, bin_frac)
data, time, exposure = datacube(directory, date, wavelength, coords, coords_type, bin_frac)


# verify each processor received subcube with correct dimensions
#ss = np.shape(subcube)  # Validation	
#print "Processor", rank, "received an array with dimensions", ss  # Validation
#print "Height = %i, Width = %i, Total pixels = %i" % (subcube.shape[0], subcube.shape[1], subcube.shape[0]*subcube.shape[1])

#params_T, M2_fit_T = spec_fit( subcube )  # Do something with the array
#params_T = spec_fit( subcube )  # Do something with the array
#newData_p = comm.gather(testr, root=0)  # Gather all the results
newData_p = comm.gather(data, root=0)  # Gather all the results
newData_t = comm.gather(time, root=0)  # Gather all the results
newData_e = comm.gather(exposure, root=0)  # Gather all the results

temp_max = np.zeros((len(newData_p)))


if rank == 0:
    for r in range(len(newData_p)):
        temp_max[r] = newData_p[r].shape[2]
    xmin = int(np.min(temp_max))
    for s in range(len(newData_p)):
        newData_p[s] = newData_p[s][:,:,0:xmin]   # ugh, the trimming should probably be off of both sides of the array
    num_f = newData_p[0].shape[0]
    y_len = newData_p[0].shape[1]
    overlap = np.zeros((num_f,20,xmin))
    full_arr = np.zeros((num_f,(len(newData_p)*y_len)-(34*(len(newData_p)-1)),xmin))
    
    for t in range(len(newData_p)-1):
        if t == 0:
            full_arr[:,0:(t+1)*(y_len-34),:] = newData_p[t][:,0:y_len-34,:]
            full_arr[:,(t+1)*(y_len-34):(t+1)*(y_len)-(t*34),:] = (newData_p[t][:,y_len-34:y_len,:] + newData_p[t+1][:,0:34,:]) / 2.
        else:
            full_arr[:,(t*y_len)-((t-1)*34):(t+1)*(y_len-34),:] = newData_p[t][:,34:y_len-34,:]
            full_arr[:,(t+1)*(y_len-34):(t+1)*(y_len)-(t*34),:] = (newData_p[t][:,y_len-34:y_len,:] + newData_p[t+1][:,0:34,:]) / 2.
    full_arr[:,(t+1)*y_len-(t*34):(t+2)*y_len-((t+1)*34), :] = newData_p[len(newData_p)-1][:,34:y_len,:] 
        
    #stack_p = np.hstack(newData_p)
    #print stack_p.shape
    print "Derotated cube shape: (%s,%s,%s)" % full_arr.shape
    #np.save('C:/Users/Brendan/Desktop/derotate_mpi_8seg_D.npy',full_arr)
    #np.save('C:/Users/Brendan/Desktop/time.npy', newData_t[0])
    #np.save('C:/Users/Brendan/Desktop/exposure.npy', newData_e[0])
    np.save('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength), full_arr)
    np.save('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength), newData_t[0])
    np.save('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength), newData_e[0])
    
    # calculate the average-intensity image of the timeseries 
    AVG = np.average(full_arr,axis=0)  #hmm, this is not normalized - doesn't really matter I don't think
    
    # determine the middle file of the timeseries
    mid_num = (full_arr.shape[0]/2)
    mid = full_arr[mid_num]
    
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
    #np.save('C:/Users/Brendan/Desktop/visual.npy', visual)
    np.save('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength), visual)

    
    T_final = timer() - start
    T_min_final, T_sec_final = divmod(T_final, 60)
    T_hr_final, T_min_final = divmod(T_min_final, 60)
    print "Total program time = %i sec" % T_final   