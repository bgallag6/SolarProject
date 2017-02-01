# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 18:38:48 2017

@author: Brendan
"""

"""
######################
# run with:
# $ mpiexec -n 8 python FFT_avg_mpi.py    (8 = number of processors)
######################
"""


"""
############################
############################
# FFT segment averaging + 3x3 Pixel Box Averaging (int)
############################
############################
"""

# added 1/29:
# if choose to save datacube as type-int, plus exposure array
# this normalizes by exposure time when extracting pixel-intensity values
# should reduce datacube array size by 4x, and not cost much time here

# not working, seems slower than usual somehow, when should be 4x quicker

# including exposure normalization in loop cause time to go nuts??

# ran through revised version, program time down from ~235 seconds to 70 seconds
# check if identical results... ugh, not exactly the same, not sure how
# off by max of 0.02, which is actually significant, not sure how happened?
# check heatmaps to see if is significant

# not sure if I did it correctly, but percent difference in some values was showing 50-70%?!?
# !! oh, possibly because I deleted two files that were odd sizes -- damn, can't compare

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
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
#import astropy.units as u
import h5py
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
import accelerate  # switch on if computer has installed
from mpi4py import MPI

#def fft_avg(subcube, timeseries, exposure_array, num_seg):
def fft_avg(subcube, timeseries, num_seg):
    
    from scipy import fftpack
    
    DATA = subcube
    
    TIME = timeseries
    
    #Ex = exposure_array
    
    #print DATA.shape 
    
    #print "Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0])
    
    t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters
    
    n_segments = num_seg  # break data into 12 segments of equal length
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
    spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))
    
    #print "length time-interp array = %i" % n
    #print "size for FFT to consider = %i" % freq_size
    #print "length of sample freq array = %i" % len(sample_freq)
    #print "length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs)
    
    for ii in range(0,spectra_seg.shape[0]):
    #for ii in range(142,160):
    
        for jj in range(0,spectra_seg.shape[1]):
        #for jj in range(0,574):        
        
            x1_box = 0+ii
            #x2_box = 2+ii  # if want to use median of more than 1x1 pixel box
            y1_box = 0+jj
            #y2_box = 2+jj  # if want to use median of more than 1x1 pixel box
            
            for k in range(0,DATA.shape[0]):
              #im=DATA[k]/(Ex[k])	  # get image + normalize by exposure time  (time went nuts?)
              im=DATA[k]	  # get image
              #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])  # finds pixel-box median
              pixmed[k]= im[x1_box,y1_box]	# median  <-- use this
            
            #pixmed = pixmed/Ex  # normalize by exposure time                
                        
            # The derotation introduces some bad data towards the end of the sequence. This trims that off
            bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
            last_good_pos = bad - 1			# retain only data before the <=zero
            
            # Get time and pixel values
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
              #sig_fft = fftpack.fft(sig)
              #sig_fft = fftpack.rfft(sig)  # real-FFT
              #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
              sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
              #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
              powers = np.abs(sig_fft)[pidxs]
              norm = len(sig)  # to normalize the power
              powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
              avg_array += powers
            
            avg_array /= n_segments  # take the average of the segments
            
            avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
                       
            spectra_seg[ii][jj] = avg_array  # construct 3D array with averaged FFTs from each pixel
        
            
    return spectra_seg
    
# load data
cube = np.load('/media/brendan/My Passport/Users/Brendan/Desktop/SolarProject/datacubes/20130530_1600_2300_2600i_2200_3000j_data_rebin2.npy')
time = np.load('/media/brendan/My Passport/Users/Brendan/Desktop/SolarProject/time_arrays/SDO_20130530_1600A_2300_2600i_2200_3000j_float_time.npy')
exposure = np.load('/media/brendan/My Passport/Users/Brendan/Desktop/SolarProject/data/20130530/193/20130530_193_2300_2600i_2200_3000j_exposure.npy')
num_seg = 6


start = timer()

comm = MPI.COMM_WORLD 	# set up comms
rank = comm.Get_rank()	# Each processor gets its own "rank"
#print "Hello World from process ", rank  # DEBUG/VALIDATE


# Rank0 is the first processor. Use that to do the main chores
if rank == 0:
  size = MPI.COMM_WORLD.Get_size()		# How many processors do we have?
  chunks = np.array_split(cube, size, axis=1)		# Split the data based on no. of processors
else:
  chunks = None  # Prepare a variable on the other nodes

# Distribute the data to all nodes, defining the root node (rank=0) as distributor
subcube = comm.scatter(chunks, root=0)
ss = np.shape(subcube)  # Validation	
print "Processor", rank, "received an array with dimensions", ss  # Validation
print "Height = %i, Width = %i, Total pixels = %i" % (subcube.shape[0], subcube.shape[1], subcube.shape[0]*subcube.shape[1])
print "Estimated time remaining... "

#spectra_seg_part = fft_avg(subcube, time, exposure, num_seg)		# Do something with the array
spectra_seg_part = fft_avg(subcube, time, num_seg)		# Do something with the array
newData_s = comm.gather(spectra_seg_part, root=0)	# Gather all the results

# Again, just have one node do the last bit
if rank == 0:
  #stack = np.vstack(newData) 	# stack the 2-d arrays together and we're done!
  spectra_seg = np.vstack(newData_s)  # should be vstack?
  print spectra_seg.shape			# Verify we have a summed version of the input cube
 
T_act = timer() - start
print "Program time = %i sec" % T_act     


# initialize arrays to hold temporary results for calculating geometric average
temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
p_geometric = np.zeros((spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
spectra_array = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
    

### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

for l in range(1,spectra_seg.shape[0]-1):
#for l in range(1,25):
    #print l
    for m in range(1,spectra_seg.shape[1]-1):
    #for m in range(1,25):
        
        temp[0] = np.log10(spectra_seg[l-1][m-1])
        temp[1] = np.log10(spectra_seg[l-1][m])
        temp[2] = np.log10(spectra_seg[l-1][m+1])
        temp[3] = np.log10(spectra_seg[l][m-1])
        temp[4] = np.log10(spectra_seg[l][m])
        temp[5] = np.log10(spectra_seg[l][m+1])
        temp[6] = np.log10(spectra_seg[l+1][m-1])
        temp[7] = np.log10(spectra_seg[l+1][m])
        temp[8] = np.log10(spectra_seg[l+1][m+1])

        temp9 = np.sum(temp, axis=0)
        p_geometric = temp9 / 9.
        spectra_array[l-1][m-1] = np.power(10,p_geometric)

np.save('/media/brendan/My Passport/Users/Brendan/Desktop/SolarProject/data/20130530/1600_2300_2600i_2200_3000j_data_rebin2_spectra_mpi', spectra_array)
