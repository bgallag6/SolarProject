# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 00:22:09 2017

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
from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate  # put inside function
import glob
import matplotlib.pyplot as plt

DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20120923/171/derotated.npy')

TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20120923/171/time.npy')

Ex = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20120923/171/exposure.npy')

t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters

#data_sub = DATA[:,149:152,400:410]

#DATA = DATA[:, 150:200, 250:300]
DATA = DATA[:, 150:151, 250:251]


"""
n_segments = 3  # break data into # segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) / n_segments


time_step = 12  # 12 second cadence 
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]
        
pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))

#for ii in range(0,DATA.shape[1]):
for ii in range(0,DATA.shape[1]):

    #for jj in range(0,DATA.shape[2]):
    for jj in range(0,DATA.shape[2]):        
    
        x1_box = 0+ii
        y1_box = 0+jj
        
        for k in range(0,DATA.shape[0]):
          im=DATA[k]
          pixmed[k]= im[x1_box,y1_box]	# median  <-- use this

        pixmed = pixmed/Ex  # normalize by exposure time    
        
        # The derotation introduces some bad data towards the end of the sequence. This trims that off
        bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
        last_good_pos = bad - 1			# retain only data before the <=zero
        
        # Get time and pixel values
        v=pixmed[0:last_good_pos]		
        t=TIME[0:last_good_pos]
        
        n_segments = 3  # break data into # segments of equal length
        n = len(t_interp)
        rem = n % n_segments
        freq_size = (n - rem) / n_segments
        
        
        time_step = 12  # 12 second cadence 
        sample_freq = fftpack.fftfreq(freq_size, d=time_step)
        pidxs = np.where(sample_freq > 0)
        freqs = sample_freq[pidxs]
               
        

        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
        data = v_interp
        
        avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers

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
        
        avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
                   
        spectra_seg[ii][jj] = avg_array  # construct 3D array with averaged FFTs from each pixel


# initialize arrays to hold temporary results for calculating arithmetic average (changed from geometric)
temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
p_avg = np.zeros((spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
spectra_array = np.zeros((spectra_seg.shape[0]-2,spectra_seg.shape[1]-2, len(freqs)))
    

### calculate 3x3 pixel-box arithmetic average.  start at 1 and end 1 before to deal with edges.
## previously for geometric -- 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

#for l in range(1,spectra_seg.shape[0]-1):
for l in range(1,DATA.shape[1]-1):
    #print l
    #for m in range(1,spectra_seg.shape[1]-1):
    for m in range(1,DATA.shape[2]-1):
                   
        temp[0] = spectra_seg[l-1][m-1]
        temp[1] = spectra_seg[l-1][m]
        temp[2] = spectra_seg[l-1][m+1]
        temp[3] = spectra_seg[l][m-1]
        temp[4] = spectra_seg[l][m]
        temp[5] = spectra_seg[l][m+1]
        temp[6] = spectra_seg[l+1][m-1]
        temp[7] = spectra_seg[l+1][m]
        temp[8] = spectra_seg[l+1][m+1]


        temp9 = np.sum(temp, axis=0)
        p_avg = temp9 / 9.
        #spectra_array[l-1][m-1] = np.power(10,p_geometric)
        spectra_array[l-1][m-1] = p_avg

        #plt.figure(figsize=(15,15))
        #plt.title('Pixel %ix, %iy | 3-Segment + 3x3 Pixel Box Averaging' % (m,l), y=1.01, fontsize=20)
        #plt.xlim(10**-4.5,10**-1)
        #plt.ylim(10**-6,10**0)
        #plt.loglog(freqs, p_avg)
        #plt.savefig('C:/Users/Brendan/Desktop/validation/3seg_3x3_pixel_%ix_%iy.jpeg' % (m,l))

np.save('C:/Users/Brendan/Desktop/validation/3seg_3x3_spectra_region.npy', spectra_array)        

"""
        
        



"""
n_segments = 1  # break data into # segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) / n_segments


time_step = 12  # 12 second cadence 
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]
        
pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))

#for ii in range(0,DATA.shape[1]):
for ii in range(0,DATA.shape[1]):

    #for jj in range(0,DATA.shape[2]):
    for jj in range(0,DATA.shape[2]):        
              
    
        x1_box = 0+ii
        y1_box = 0+jj
        
        for k in range(0,DATA.shape[0]):
          im=DATA[k]
          pixmed[k]= im[x1_box,y1_box]	# median  <-- use this

        pixmed = pixmed/Ex  # normalize by exposure time    
        
        # The derotation introduces some bad data towards the end of the sequence. This trims that off
        bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
        last_good_pos = bad - 1			# retain only data before the <=zero
        
        # Get time and pixel values
        v=pixmed[0:last_good_pos]		
        t=TIME[0:last_good_pos]
        
        n_segments = 1  # break data into # segments of equal length
        n = len(t_interp)
        rem = n % n_segments
        freq_size = (n - rem) / n_segments
        
        
        time_step = 12  # 12 second cadence 
        sample_freq = fftpack.fftfreq(freq_size, d=time_step)
        pidxs = np.where(sample_freq > 0)
        freqs = sample_freq[pidxs]
        
    
        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
        data = v_interp
                   
        ## perform Fast Fourier Transform on each segment       
        sig = data
        sig_fft = fftpack.fft(sig)
        #sig_fft = fftpack.rfft(sig)  # real-FFT
        #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
        #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
        #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
        powers = np.abs(sig_fft)[pidxs]
        norm = len(sig)  # to normalize the power
        powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
        
        avg_array = np.transpose(powers)  # take transpose of array to fit more cleanly in 3D array
                   
        spectra_seg[ii][jj] = avg_array  # construct 3D array with averaged FFTs from each pixel


# initialize arrays to hold temporary results for calculating arithmetic average (changed from geometric)
temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
p_avg = np.zeros((spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
spectra_array = np.zeros((spectra_seg.shape[0]-2,spectra_seg.shape[1]-2, len(freqs)))
    

### calculate 3x3 pixel-box arithmetic average.  start at 1 and end 1 before to deal with edges.
## previously for geometric -- 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

#for l in range(1,spectra_seg.shape[0]-1):
for l in range(1,DATA.shape[1]-1):
    #print l
    #for m in range(1,spectra_seg.shape[1]-1):
    for m in range(1,DATA.shape[2]-1):
                   
        temp[0] = spectra_seg[l-1][m-1]
        temp[1] = spectra_seg[l-1][m]
        temp[2] = spectra_seg[l-1][m+1]
        temp[3] = spectra_seg[l][m-1]
        temp[4] = spectra_seg[l][m]
        temp[5] = spectra_seg[l][m+1]
        temp[6] = spectra_seg[l+1][m-1]
        temp[7] = spectra_seg[l+1][m]
        temp[8] = spectra_seg[l+1][m+1]


        temp9 = np.sum(temp, axis=0)
        p_avg = temp9 / 9.
        #spectra_array[l-1][m-1] = np.power(10,p_geometric)
        spectra_array[l-1][m-1] = p_avg

        #plt.figure(figsize=(15,15))
        #plt.title('Pixel %ix, %iy | 1-Segment + 3x3 Pixel Box Averaging' % (m,l), y=1.01, fontsize=20)
        #plt.xlim(10**-5,10**-1)
        #plt.ylim(10**-6,10**0)
        #plt.loglog(freqs, p_avg)
        #plt.savefig('C:/Users/Brendan/Desktop/validation/1seg_3x3_pixel_%ix_%iy.jpeg' % (m,l))
        
np.save('C:/Users/Brendan/Desktop/validation/1seg_3x3_spectra_region.npy', spectra_array)
"""




#"""
n_segments = 1  # break data into # segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) / n_segments


time_step = 12  # 12 second cadence 
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]
        
pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))

spectra_array = np.zeros((spectra_seg.shape[0],spectra_seg.shape[1], len(freqs)))

#for ii in range(0,DATA.shape[1]):
for ii in range(0,DATA.shape[1]):

    #for jj in range(0,DATA.shape[2]):
    for jj in range(0,DATA.shape[2]):        
             
    
        x1_box = 0+ii
        y1_box = 0+jj
        
        for k in range(0,DATA.shape[0]):
          im=DATA[k]
          pixmed[k]= im[x1_box,y1_box]	# median  <-- use this

        pixmed = pixmed/Ex  # normalize by exposure time    
        
        # The derotation introduces some bad data towards the end of the sequence. This trims that off
        bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
        last_good_pos = bad - 1			# retain only data before the <=zero
        
        # Get time and pixel values
        v=pixmed[0:last_good_pos]		
        t=TIME[0:last_good_pos]
        
        n_segments = 1  # break data into # segments of equal length
        n = len(t_interp)
        rem = n % n_segments
        freq_size = (n - rem) / n_segments
        
        
        time_step = 12  # 12 second cadence 
        sample_freq = fftpack.fftfreq(freq_size, d=time_step)
        pidxs = np.where(sample_freq > 0)
        freqs = sample_freq[pidxs]
        
    
        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
        data = v_interp
        
        #norm_v = (v_interp-np.min(v_interp)) / (np.max(v_interp) - np.min(v_interp)) 
        avg_v = np.average(v_interp)
        norm_avg = (v_interp - avg_v) / avg_v
        
        plt.figure(figsize=(15,15))
        plt.title('1-Segment + No Averaging -- Timeseries Normalized', y=1.01, fontsize=20)
        plt.plot(t_interp,norm_avg)
        plt.savefig('C:/Users/Brendan/Desktop/validation/1seg_1x1_timeseries_normalized.jpeg')
        
        #plt.figure()
        #plt.plot(t_interp,norm_avg)
                   
        ## perform Fast Fourier Transform on each segment       
        sig = norm_avg
        sig_fft = fftpack.fft(sig)
        #sig_fft = fftpack.rfft(sig)  # real-FFT
        #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
        #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
        #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
        powers = np.abs(sig_fft)[pidxs]
        norm = len(sig)  # to normalize the power
        powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
        
        spectra_array[ii][jj] = powers


        plt.figure(figsize=(15,15))
        plt.title('1-Segment + No Averaging -- Timeseries Normalized', y=1.01, fontsize=20)
        plt.xlim(10**-5,10**-1)
        plt.ylim(10**-6,10**0)
        plt.loglog(freqs, powers)
        plt.savefig('C:/Users/Brendan/Desktop/validation/1seg_1x1_timeseries_normalized_spectra.jpeg')

#np.save('C:/Users/Brendan/Desktop/validation/1seg_1x1_spectra_region.npy', spectra_array)        

#"""
        
        
        
        
"""
n_segments = 3  # break data into # segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) / n_segments


time_step = 12  # 12 second cadence 
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]
        
pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))

spectra_array = np.zeros((spectra_seg.shape[0],spectra_seg.shape[1], len(freqs)))

#for ii in range(0,DATA.shape[1]):
for ii in range(0,DATA.shape[1]):

    #for jj in range(0,DATA.shape[2]):
    for jj in range(0,DATA.shape[2]):        
             
    
        x1_box = 0+ii
        y1_box = 0+jj
        
        for k in range(0,DATA.shape[0]):
          im=DATA[k]
          pixmed[k]= im[x1_box,y1_box]	# median  <-- use this

        pixmed = pixmed/Ex  # normalize by exposure time    
        
        # The derotation introduces some bad data towards the end of the sequence. This trims that off
        bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
        last_good_pos = bad - 1			# retain only data before the <=zero
        
        # Get time and pixel values
        v=pixmed[0:last_good_pos]		
        t=TIME[0:last_good_pos]
        
        n_segments = 3  # break data into # segments of equal length
        n = len(t_interp)
        rem = n % n_segments
        freq_size = (n - rem) / n_segments
        
        
        time_step = 12  # 12 second cadence 
        sample_freq = fftpack.fftfreq(freq_size, d=time_step)
        pidxs = np.where(sample_freq > 0)
        freqs = sample_freq[pidxs]
        
    
        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
        data = v_interp
        
        avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers

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
        
        avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array

        spectra_array[ii][jj] = avg_array 
        
        #plt.figure(figsize=(15,15))
        #plt.title('Pixel %ix, %iy | 3-Segment + 1x1 Pixel Box Averaging' % (jj,ii), y=1.01, fontsize=20)
        #plt.xlim(10**-5,10**-1)
        #plt.ylim(10**-6,10**0)
        #plt.loglog(freqs, avg_array)
        #plt.savefig('C:/Users/Brendan/Desktop/validation/3seg_1x1_pixel_%ix_%iy.jpeg' % (jj,ii))

np.save('C:/Users/Brendan/Desktop/validation/3seg_1x1_spectra_region.npy', spectra_array)       
"""