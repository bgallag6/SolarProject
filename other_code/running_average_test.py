# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 16:04:24 2017

@author: Brendan
"""
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Cursor
from pylab import *
import numpy as np
import sunpy
from sunpy.map import Map
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate  # switch on if computer has installed




from scipy import fftpack

DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/193/20130626_193_-450_-200i_-200_200j_data_rebin1.npy')
TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/193/20130626_193_-450_-200i_-200_200j_time.npy')
EXPOSURE = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/193/20130626_193_-450_-200i_-200_200j_exposure.npy')


Ex = EXPOSURE

print DATA.shape 

print "Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0])

t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters


num_seg = 1

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

print "length time-interp array = %i" % n
print "size for FFT to consider = %i" % freq_size
print "length of sample freq array = %i" % len(sample_freq)
print "length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs)


start = timer()
T1 = 0

#for ii in range(0,spectra_seg.shape[0]):
for ii in range(100,101):

    #for jj in range(0,spectra_seg.shape[1]):
    for jj in range(100,101):        
    
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
        v=pixmed[0:last_good_pos]		
        t=TIME[0:last_good_pos]
        #v=pixmed  # use for 335/131/094 -- can't get rid of negative values for those
        #t=TIME
        
        #plt.plot(t,v)
    
        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
        v_interp_avg = np.zeros_like(v_interp)
        v_interp_average = 0          
        
        for w in range(len(v_interp)):
            if w == 0:
                v_interp_avg[w] = v_interp[w]
            elif w == 1:
                v_interp_avg[w] = v_interp[w]
            elif w == len(v_interp)-2:
                v_interp_avg[w] = v_interp[w]
            elif w == len(v_interp)-1:
                v_interp_avg[w] = v_interp[w]
            else:
                v_interp_average = v_interp[w-2]+v_interp[w-1]+v_interp[w]+v_interp[w+1]+v_interp[w+2]
                v_interp_avg[w] = v_interp_average / 5.
        
        
        data = v_interp_avg
        plt.figure()
        plt.plot(t_interp,v_interp)
        plt.figure()
        plt.plot(t_interp,v_interp_avg)
        
        avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
        
                
        sig = data            
        sig_fft = fftpack.fft(sig)
        powers = np.abs(sig_fft)[pidxs]
        norm = len(sig)  # to normalize the power
        powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
        
        plt.figure()
        plt.xlim(10**-5,10**-1)
        plt.ylim(10**-8,10**0)
        plt.loglog(freqs, powers)
        
        
        data = v_interp
        data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
        split = np.split(data, n_segments)  # create split array for each segment        
       
        for i in range(0,n_segments):               
            
          ## perform Fast Fourier Transform on each segment       
          sig = split[i]
          sig_fft = fftpack.fft(sig)
          #sig_fft = fftpack.rfft(sig)  # real-FFT
          #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
          #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
          #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
          powers = np.abs(sig_fft)[pidxs]
          norm = len(sig)  # to normalize the power
          powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
          
       
        plt.figure()
        plt.xlim(10**-5,10**-1)
        plt.ylim(10**-8,10**0)
        plt.loglog(freqs, powers)
        
        
        
        num_seg = 6

        n_segments = num_seg  # break data into 12 segments of equal length
        n = len(t_interp)
        rem = n % n_segments
        freq_size = (n - rem) / n_segments
        
        ## determine frequency values that FFT will evaluate
        time_step = 12  # add as argument in function call, or leave in as constant?
        sample_freq = fftpack.fftfreq(freq_size, d=time_step)
        pidxs = np.where(sample_freq > 0)
        freqs = sample_freq[pidxs]
        
        avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
        
        data = v_interp
        data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
        split = np.split(data, n_segments)  # create split array for each segment        
       
        for i in range(0,n_segments):               
            
          ## perform Fast Fourier Transform on each segment       
          sig = split[i]
          sig_fft = fftpack.fft(sig)
          #sig_fft = fftpack.rfft(sig)  # real-FFT
          #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
          #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
          #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
          powers = np.abs(sig_fft)[pidxs]
          norm = len(sig)  # to normalize the power
          powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
          avg_array += powers          
        
        avg_array /= n_segments  # take the average of the segments
       
        plt.figure()
        plt.xlim(10**-5,10**-1)
        plt.ylim(10**-8,10**0)
        plt.loglog(freqs, avg_array)