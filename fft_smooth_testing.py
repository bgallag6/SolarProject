# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:44:24 2017

@author: Brendan
"""

import numpy as np
import scipy.signal
from pylab import *
from sunpy.map import Map
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate  # put inside function
import glob

import accelerate  # switch on if computer has installed

from scipy import fftpack



DATA = np.load('C:/Users/Brendan/Desktop/testing_171/derotated.npy')

TIME = np.load('C:/Users/Brendan/Desktop/testing_171/time.npy')

Ex = np.load('C:/Users/Brendan/Desktop/testing_171/exposure.npy')

wavelength = 171

print DATA.shape 

print "Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0])

t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters

n_segments = 1  # break data into 12 segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) / n_segments

## determine frequency values that FFT will evaluate
if wavelength == 1600 or wavelength == 1700:
  time_step = 24  # 24 second cadence for these wavelengths
else:
  time_step = 12  # 12 second cadence for the others
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

reslt = (DATA.shape[0] == TIME.shape[0])
print "DATA and TIME array sizes match: %s" % reslt

pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))
spectra_array = np.zeros((spectra_seg.shape[0], spectra_seg.shape[1], spectra_seg.shape[2]))  # would pre-allocating help? (seems to)

print "length time-interp array = %i" % n
print "size for FFT to consider = %i" % freq_size
print "length of sample freq array = %i" % len(sample_freq)
print "length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs)


start = timer()
T1 = 0

for ii in range(0,spectra_seg.shape[0]):
#for ii in range(0,5):

    for jj in range(0,spectra_seg.shape[1]):
    #for jj in range(0,5):        
    
        x1_box = 0+ii
        #x2_box = 2+ii  # if want to use median of more than 1x1 pixel box
        y1_box = 0+jj
        #y2_box = 2+jj  # if want to use median of more than 1x1 pixel box
        
        for k in range(0,DATA.shape[0]):
          #im=DATA[k]/(Ex[k])	  # get image + normalize by exposure time  (time went nuts?)
          im=DATA[k]
          #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])  # finds pixel-box median
          pixmed[k]= im[x1_box,y1_box]	# median  <-- use this

        pixmed = pixmed/Ex  # normalize by exposure time    
        
        # The derotation introduces some bad data towards the end of the sequence. This trims that off
        bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
        last_good_pos = bad - 1			# retain only data before the <=zero
        
        # Get time and pixel values
        if wavelength == 94 or wavelength == 131 or wavelength == 335:  
            v=pixmed  # 094/131/335 -- intensities too low to trim off negative values
            t=TIME
        else:
            v=pixmed[0:last_good_pos]		
            t=TIME[0:last_good_pos]
        
    
        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
        data = v_interp
              
        ## perform Fast Fourier Transform on each segment       
        sig = data
        #sig_fft = fftpack.fft(sig)
        #sig_fft = fftpack.rfft(sig)  # real-FFT
        #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
        sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
        #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
        powers = np.abs(sig_fft)[pidxs]
        norm = len(sig)  # to normalize the power
        powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
       
        orig_powerspec = powers

        # Now Smooth the Power Spectra
        orig = orig_powerspec.copy()		# make a copy
        
        window_len=11   # Keep this an odd number
        
        window='hanning'  # OTHER CHOICES =  'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        s=np.r_[orig[window_len-1:0:-1],orig,orig[-1:-window_len:-1]]
        
        if window == 'flat': #moving average
        	w=np.ones(window_len,'d')
        else:
        	w=eval('np.'+window+'(window_len)')
        
        	
        y=np.convolve(w/w.sum(),s,mode='valid')
        y=y[( window_len/2+1): len(orig)+( window_len/2 +1)]	# Crop down the new TS
        
        # Now smooth the smooth...
        y2c = y.copy()
        s2=np.r_[y2c[window_len-1:0:-1],y2c,y2c[-1:-window_len:-1]]
        y2=np.convolve(w/w.sum(),s2,mode='valid')
        y2=y2[( window_len/2+1): len(orig)+( window_len/2 +1)]   # Crop down the new TS
                   
        spectra_array[ii][jj] = powers  # construct 3D array with averaged FFTs from each pixel
        
    # estimate time remaining and print to screen
    T = timer()
    T2 = T - T1
    if ii == 0:
        T_init = T - start
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
T_act = timer() - start
T_min3, T_sec3 = divmod(T_act, 60)
T_hr3, T_min3 = divmod(T_min3, 60)
print "Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3) 

       
        
np.save('C:/Users/Brendan/Desktop/testing_171/spectra_1seg_raw.npy', spectra_array)


