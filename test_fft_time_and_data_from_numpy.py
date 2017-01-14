# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 13:19:01 2016

@author: Brendan
"""
"""
# need to put this into full module automated
# doesn't like when time is not an integer, get error:
# t_interp = [(12*l) for l in range(0,TIME[len(TIME)-1]/12)]
# TypeError: range() integer end argument expected, got numpy.float64.

# update 1/10 : taking out of loop - to calculate # frequencies 

# generalize so that t_interp is valid even if not starting at 00:00:00

# maybe start thinking about error checking

# took out fitting functions - obviously dont need them

# took out a bunch of module imports - dont need for this function

# full passthrough from this to 3x3 worked perfectly for one run of 1600 rebin4

# update 1/12: added estimated time remaining

# - insane, I did profiling and 90% of time was spent on finding pixel median values
# - however, it was the median of one pixel.  I took that out and the results are identical
# - in 10% of the time

# t_interp issue -- ran through full program and was exactly same 
# changed t_interp to linspace with one extra point, took out conversion to float in loop

# possible function arguments: DATA, TIME, n_segments, time_step, median-size?, custom range within subregion?

# included in module 12:30 AM 1/13
"""

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
from scipy import fftpack
from timeit import default_timer as timer



    
#with h5py.File('C:/Users/Brendan/Desktop/SDO/SDO_20120923_211A_(528)_(132)x_(100)_100y.hdf5', 'r', driver='core') as f:
    #DATA = f['Image Data']
    #TIME = f['Time Data']
    
DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/datacubes/20140902_193_2000_2900i_1850_3050j_data_rebin2.npy')

TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/time_arrays/20140902_193_2000_2900i_1850_3050j_time.npy')

print DATA.shape 

print "Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0])

t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters

n_segments = 6  # break data into 12 segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) / n_segments

## determine frequency values that FFT will evaluate
time_step = 12  # add as argument in function call, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

reslt = (DATA.shape[0] == TIME.shape[0])
print "DATA and TIME array sizes match: %s" % reslt

pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_array = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))

print "length time-interp array = %i" % n
print "size for FFT to consider = %i" % freq_size
print "length of sample freq array = %i" % len(sample_freq)
print "length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs)


start = timer()
T1 = 0

for ii in range(0,spectra_array.shape[0]):
#for ii in range(142,160):

    for jj in range(0,spectra_array.shape[1]):
    #for jj in range(0,574):        
    
        x1_box = 0+ii
        x2_box = 2+ii
        y1_box = 0+jj
        y2_box = 2+jj
        
        for k in range(0,DATA.shape[0]):
                    im=DATA[k]												# get image
                    #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])	# median
                    pixmed[k]= im[x1_box,y1_box]	# median  <-- use this
                    
                    
        # The derotation introduces some bad data towards the end of the sequence. This trims that off
        bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
        last_good_pos = bad - 1			# retain only data before the <=zero
        
        # Get t and v data
        v=pixmed[0:last_good_pos]		
        t=TIME[0:last_good_pos]
        
        #plt.plot(t,v)
    
        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
        data = v_interp
        
        avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers

        data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
        split = np.split(data, n_segments)  # create split array for each segment


        for i in range(0,n_segments):               
            
             ## perform Fast Fourier Transform on each segment       
             sig = split[i]
             sig_fft = fftpack.fft(sig)
             powers = np.abs(sig_fft)[pidxs]
             norm = len(sig)  # to normalize the power
             powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
             avg_array += powers
        
        avg_array /= n_segments  # take the average of the segments
        
        avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
       
        spectra_array[ii][jj] = avg_array  # construct 3D array with averaged FFTs from each pixel
    
    # estimate time remaining and print to screen     
    T = timer() - start
    T2 = T - start - T1
    if ii == 0:
        T_est = T2*(spectra_array.shape[0])    
    T_est2 = T2*(spectra_array.shape[0]-ii)
    print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_array.shape[0], T_est2)
    T1 = T

# print estimated and total program time to screen 
print "Beginning Estimated time = %i sec" % T_est
T_act = timer() - start
print "Actual total time = %i sec" % T_act 