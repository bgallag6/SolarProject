# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 18:38:30 2017

@author: Brendan
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pylab import *
import glob
import sunpy
from sunpy.map import Map
#from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
import numpy as np
import astropy.units as u
from sunpy.image.coalignment import calculate_match_template_shift
from sunpy.physics.transforms.solar_rotation import calculate_solar_rotate_shift


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
# or possibly import these from the return values of arc2pix  <-- this
x1 = -500
x2 = 500
y1 = -400
y2 = 400

bin_frac = 1

# create a list of all the files. This is USER-DEFINED
flist = glob.glob('F:/Users/Brendan/Desktop/SolarProject/data/20130530/1712/aia*.fits')
nf = len(flist)

# Select the image that is the "middle" of our selection.
# We do this because the solar derotation algorithm operates centered on the 
# "middle" image  (might not be used)
mid_file = np.int(np.floor(nf / 2))

 

mc_list = []  # create an empty list

count = 0  # counter to see where program is at 

# Use defined coordinates, extract the submaps from each AIA image, and store
# them in the empty list. This takes many minutes to complete.
print " "
print "Reading files and extracting submaps. This takes a while..."
print " "
   



for filename in flist:
    mc_list.append(Map(filename).submap([x1,x2]*u.arcsecond, [y1,y2]*u.arcsecond))
    if (count%10) == 0:
       print("file %i out of %i" % (count,nf))  # just to see where program is at
    count += 1
    
	
# Create a new Map Cube object to hold our de-rotated data
new_mapcube = Map(mc_list, cube=True)
#ani = new_mapcube.plot()   # doctest: +SKIP
#ani = new_mapcube.peek()   # doctest: +SKIP
#plt.show()   # doctest: +SKIP

shifts = calculate_solar_rotate_shift(new_mapcube, layer_index=mid_file)
#print shifts
L = shifts['x']
L = np.array(L)

# Perform the derotation of the submaps. This take a while too.
dr = mapcube_solar_derotate(new_mapcube, layer_index=0)
#dr = mapcube_solar_derotate(new_mapcube, layer_index=6)
#dr = mapcube_solar_derotate(new_mapcube, layer_index=11)

#ani = dr.plot()
#plt.show()

print "done derotating"


mid_subarr = dr[mid_file].data		# extract data from middle file of derotated datacube
rem_i = mid_subarr.shape[0] % bin_frac  # calculate remainder to get integer dimensions
rem_j = mid_subarr.shape[1] % bin_frac  # calculate remainder to get integer dimensions
subarr_idim = (mid_subarr.shape[0] - rem_i) / bin_frac  # get ydim
subarr_jdim = (mid_subarr.shape[1] - rem_j) / bin_frac  # get xdim
 

t = dr[0].date  # extract the date / time from the first image
base_time = (t.hour * 60.) + (t.minute * 60.) + t.second  # convert date / time to seconds

# initialize arrays to hold exposure time, pixel data, and time values
I = np.empty((nf))  # exposure time
DATA = np.empty((nf, subarr_idim, subarr_jdim), dtype=np.int16)  # save as int16, since that is what original is
#TIME = np.empty((nf), dtype='int')
TIME = np.empty((nf))  # might as well just have all as float

# loop through datacube and extract pixel data and time values
for p in range(0,nf):
    Ex = dr[p].exposure_time
    I[p] = Ex.value
    L = dr[p].data
    L_trim = L[0:(mid_subarr.shape[0] - rem_i), 0:(mid_subarr.shape[1] - rem_j)]
    small_L = rebin(L_trim, L_trim.shape[0]/bin_frac, L_trim.shape[1]/bin_frac)
    #DATA[p][:][:] = small_L/Ex  # normalize by exposure time
    DATA[p][:][:] = small_L  # normalize by exposure time
    T = dr[p].date
    curr_time=(T.hour * 3600.)+(T.minute * 60.)+T.second	
    TIME[p] = curr_time - base_time  # calculate running time of image

"""
for i in range(12):
    fig = plt.figure()
    plt.imshow(DATA[i])
"""
# calculate the average-intensity image of the timeseries 
#AVG = np.average(DATA,axis=0)


#fig = plt.figure(figsize=(20,20))
#plt.imshow(AVG)



#shifts = calculate_match_template_shift(new_mapcube, layer_index=1)
shifts = calculate_solar_rotate_shift(new_mapcube, layer_index=mid_file)
#print shifts
L = shifts['x']
L = np.array(L)
#M = shifts['y']
#M = np.array(M)
print L
#print M

"""
x_coord = m1.data_to_pixel(x1*u.arcsec,x2*u.arcsec,1)
x1_coord = x_coord[0].value
x2_coord = x_coord[1].value
y_coord = m1.data_to_pixel(y1*u.arcsec,y2*u.arcsec,1)
y1_coord = y_coord[0].value
y2_coord = y_coord[1].value
"""
#print data_to_pixel()