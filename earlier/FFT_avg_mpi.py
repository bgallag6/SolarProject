# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 18:38:48 2017

@author: Brendan
"""

"""
######################
# run with:
# $ mpiexec -n # python FFT_avg_mpi.py    (# = number of processors)
######################
"""


"""
############################
############################
# FFT segment averaging + 3x3 Pixel Box Averaging
############################
############################
"""


import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *
from scipy import signal
import scipy.misc
from timeit import default_timer as timer
#import accelerate  # switch on if computer has installed
from mpi4py import MPI
from scipy import fftpack

def fft_avg(subcube):
    
    from scipy import fftpack
    
    DATA = subcube
    
    reslt = (DATA.shape[0] == TIME.shape[0])
    print("DATA and TIME array sizes match: %s" % reslt, flush=True)
    
    pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
    spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))
    
    print("length time-interp array = %i" % n, flush=True)
    print("length of sample freq array = %i" % len(sample_freq), flush=True)
    print("length of freqs array = %i (should be 1/2 of row above)" % len(freqs), flush=True)
    
    start_sub = timer()
    T1 = 0    
    
    for ii in range(0,spectra_seg.shape[0]):
    
        for jj in range(0,spectra_seg.shape[1]):        
        
        
            pixmed = DATA[:,ii,jj] / exposure  # extract timeseries + normalize by exposure time   
        
            v_interp = np.interp(t_interp,TIME,pixmed)  # interpolate pixel-intensity values onto specified time grid
            
            data = v_interp
            
            avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
    
            data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
            split = np.split(data, n_segments)  # create split array for each segment
       
            for i in range(0,n_segments):               
                
              ## perform Fast Fourier Transform on each segment       
              sig = split[i]
              sig_fft = fftpack.fft(sig)
              #sig_fft = fftpack.rfft(sig)  # real-FFT                
              #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
              #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
              powers = np.abs(sig_fft)[pidxs]
              norm = len(sig)
              powers = ((powers/norm)**2)*(1./(sig.std()**2))*2   # normalize the power
              avg_array += powers
            
            avg_array /= n_segments  # take the average of the segments
            
            avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
                       
            spectra_seg[ii][jj] = avg_array  # construct 3D array with averaged FFTs from each pixel
        

            
        # estimate time remaining and print to screen
        T = timer()
        T2 = T - T1
        if ii == 0:
            T_init = T - start_sub
            T_est = T_init*(spectra_seg.shape[0])  
            T_min, T_sec = divmod(T_est, 60)
            T_hr, T_min = divmod(T_min, 60)
            #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est)
            print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spectra_seg.shape[0], T_hr, T_min, T_sec), flush=True)
        else:
            T_est2 = T2*(spectra_seg.shape[0]-ii)
            T_min2, T_sec2 = divmod(T_est2, 60)
            T_hr2, T_min2 = divmod(T_min2, 60)
            #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est2)
            print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spectra_seg.shape[0], T_hr2, T_min2, T_sec2), flush=True)
        T1 = T
        
        
    # print estimated and total program time to screen 
    print("Beginning Estimated time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec), flush=True)
    T_act = timer() - start_sub
    T_min3, T_sec3 = divmod(T_act, 60)
    T_hr3, T_min3 = divmod(T_min3, 60)
    print("Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3), flush=True) 
                    
    return spectra_seg


"""
# Setup MPI and load datacube, time, and exposure arrays
"""  
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
	
start = timer()

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command) 


import sys

# set variables from command line
directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])
  
cube = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))
#cube = np.memmap('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_data_rebin1_mmap.npy', dtype='int16', mode='r', shape=(2926,297,630))

chunks = np.array_split(cube, size, axis=1)  # Split the data based on no. of processors

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]

TIME = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))
exposure = np.load('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength))
num_seg = 6

# determine frequency values that FFT will evaluate
if wavelength in [1600,1700]:
    time_step = 24  # add as argument in function call, or leave in as constant?
else:
    time_step = 12

t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/time_step)+1)  # interpolate onto default-cadence time-grid
    
n_segments = num_seg  # break data into # segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) // n_segments

sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]



"""
### Each processor runs function on subcube, results are gathered when finished
"""
# verify each processor received subcube with correct dimensions
ss = np.shape(subcube)  # Validation	
print("Processor", rank, "received an array with dimensions", ss, flush=True)  # Validation

spectra_seg_part = fft_avg(subcube)  # Do something with the array
newData_s = comm.gather(spectra_seg_part, root=0)  # Gather all the results

# Have one node do the last bit
if rank == 0:
    spectra_seg = np.vstack(newData_s)
    print(spectra_seg.shape, flush=True)  # Verify we have a summed version of the input cube
 

    """
    ### 3x3 Averaging
    """

    temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
    spectra_array = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))
    spectra_StdDev = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))
    
    ### calculate 3x3 pixel-box average, start at 1 and end 1 before to deal with edges
    
    for l in range(1,spectra_seg.shape[0]-1):
    
        for m in range(1,spectra_seg.shape[1]-1):
            
            temp[0] = spectra_seg[l-1][m-1]
            temp[1] = spectra_seg[l-1][m]
            temp[2] = spectra_seg[l-1][m+1]
            temp[3] = spectra_seg[l][m-1]
            temp[4] = spectra_seg[l][m]
            temp[5] = spectra_seg[l][m+1]
            temp[6] = spectra_seg[l+1][m-1]
            temp[7] = spectra_seg[l+1][m]
            temp[8] = spectra_seg[l+1][m+1]
    
            p_avg = np.average(temp, axis=0)
            
            spectra_array[l-1][m-1] = p_avg
            spectra_StdDev[l-1][m-1] = np.std(temp, axis=0)
    
    T_final = timer() - start
    T_min_final, T_sec_final = divmod(T_final, 60)
    T_hr_final, T_min_final = divmod(T_min_final, 60)
    print("Total program time = %i sec" % T_final, flush=True)
    

    np.save('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength), spectra_array)
    np.save('%s/DATA/Temp/%s/%i/3x3_stddev.npy' % (directory, date, wavelength), spectra_StdDev)
