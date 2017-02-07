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

# including exposure normalization in loop cause time to go nuts??

# ran through revised version, program time down from ~235 seconds to 70 seconds
# check if identical results... ugh, not exactly the same, not sure how
# off by max of 0.02, which is actually significant, not sure how happened?
# check heatmaps to see if is significant

# not sure if I did it correctly, but percent difference in some values was showing 50-70%?!?
# !! oh, possibly because I deleted two files that were odd sizes -- damn, can't compare

# update 2/1:
# took time, exposure, freqs out of loop, since same for all.

# time for 300x600 20130530 193 region w/accelerate = 82 seconds 
# scaled to 1500x1500, should take max of 1025 seconds = 17 minutes w/ 4 cores

# for file naming / calling : put all in folder - glob fname list.  
# if certain year, certain wavelength, if last letters are spectra, pull that...(handle non-exact region sizes)
# specify path, date, wavelength - thats it.  

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
from timeit import default_timer as timer
import accelerate  # switch on if computer has installed
from mpi4py import MPI
from scipy import fftpack

def fft_avg(subcube):
    
    from scipy import fftpack
    
    DATA = subcube
    
    reslt = (DATA.shape[0] == TIME.shape[0])
    print "DATA and TIME array sizes match: %s" % reslt
    
    pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
    spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))
    
    print "length time-interp array = %i" % n
    print "length of sample freq array = %i" % len(sample_freq)
    print "length of freqs array = %i (should be 1/2 of row above)" % len(freqs)
    
    start_sub = timer()
    T1 = 0    
    
    for ii in range(0,spectra_seg.shape[0]):
    #for ii in range(0,50):
    
        for jj in range(0,spectra_seg.shape[1]):
        #for jj in range(0,200):        
        
            x1_box = 0+ii
            #x2_box = 2+ii  # if want to use median of more than 1x1 pixel box
            y1_box = 0+jj
            #y2_box = 2+jj  # if want to use median of more than 1x1 pixel box
            
            for k in range(0,DATA.shape[0]):
              #im=DATA[k]/(Ex[k])	  # get image + normalize by exposure time  (time went nuts?)
              im=DATA[k]	  # get image
              #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])  # finds pixel-box median
              pixmed[k]= im[x1_box,y1_box]	# median  <-- use this
            
            pixmed = pixmed/exposure  # normalize by exposure time                
                        
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
            
        # estimate time remaining and print to screen
        T = timer()
        T2 = T - T1
        if ii == 0:
            T_init = T - start_sub
            T_est = T_init*(spectra_seg.shape[0])  
            T_min, T_sec = divmod(T_est, 60)
            T_hr, T_min = divmod(T_min, 60)
            #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est)
            print "Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spectra_seg.shape[0], T_hr, T_min, T_sec)
        else:
            T_est2 = T2*(spectra_seg.shape[0]-ii)
            T_min2, T_sec2 = divmod(T_est2, 60)
            T_hr2, T_min2 = divmod(T_min2, 60)
            #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est2)
            print "Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spectra_seg.shape[0], T_hr2, T_min2, T_sec2)
        T1 = T
        
        
    # print estimated and total program time to screen 
    print "Beginning Estimated time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec)
    T_act = timer() - start_sub
    T_min3, T_sec3 = divmod(T_act, 60)
    T_hr3, T_min3 = divmod(T_min3, 60)
    print "Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3) 
                    
    return spectra_seg


"""
# Setup MPI and load datacube, time, and exposure arrays
"""  
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
	
start = timer()

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command) 
  
cube = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_data_rebin1.npy')
#cube = np.memmap('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_data_rebin1_mmap.npy', dtype='int16', mode='r', shape=(2926,297,630))

chunks = np.array_split(cube, size, axis=1)  # Split the data based on no. of processors

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]

TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_time.npy')
exposure = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_exposure.npy')
num_seg = 6

t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  # interpolate time array onto uniform 12-second grid
    
n_segments = num_seg  # break data into # segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) / n_segments

# determine frequency values that FFT will evaluate
time_step = 12  # add as argument in function call, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]


"""
### Each processor runs function on subcube, results are gathered when finished
"""
# verify each processor received subcube with correct dimensions
ss = np.shape(subcube)  # Validation	
print "Processor", rank, "received an array with dimensions", ss  # Validation
print "Height = %i, Width = %i, Total pixels = %i" % (subcube.shape[1], subcube.shape[2], subcube.shape[1]*subcube.shape[2])

spectra_seg_part = fft_avg(subcube)  # Do something with the array
newData_s = comm.gather(spectra_seg_part, root=0)  # Gather all the results

# Again, just have one node do the last bit
if rank == 0:
  spectra_seg = np.vstack(newData_s)  # should be vstack?
  print spectra_seg.shape  # Verify we have a summed version of the input cube
 


"""
### 3x3 Averaging
"""
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

T_final = timer() - start
T_min_final, T_sec_final = divmod(T_final, 60)
T_hr_final, T_min_final = divmod(T_min_final, 60)
print "Total program time = %i sec" % T_final

#np.save('C:/Users/Brendan/Desktop/20130530_193_2300_2600i_2200_3000j_spectra_mpi_final', spectra_array)
