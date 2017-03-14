# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 14:46:21 2017

@author: Brendan
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
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate  # switch on if computer has installed



from scipy import fftpack

DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193/derotated.npy')
TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193/time.npy')
Ex = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193/exposure.npy')

num_seg = 6

print DATA.shape 

print "Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0])

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

print "length time-interp array = %i" % n
print "size for FFT to consider = %i" % freq_size
print "length of sample freq array = %i" % len(sample_freq)
print "length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs)


start = timer()
T1 = 0

#for ii in range(0,spectra_seg.shape[0]):
for ii in range(11,12):

    #for jj in range(0,spectra_seg.shape[1]):
    for jj in range(11,12):        
    
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
        
        #data = v_interp
        
        avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers

        #data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
        #data2 = data[(len(data)/n_segments)/2:(len(data)-((len(data)/n_segments)/2))]
        #split = np.split(data, n_segments)  # create split array for each segment
        #split2 = np.split(data2, n_segments-1)

        """
        for i in range(0,n_segments):               
          if i < (n_segments-1): 
              ## perform Fast Fourier Transform on each segment       
              sig = split[i]
              #sig2 = split2[i]
              sig_fft = fftpack.fft(sig)
              #sig2_fft = fftpack.fft(sig2)
              #sig_fft = fftpack.rfft(sig)  # real-FFT
              #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
              #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
              #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
              powers = np.abs(sig_fft)[pidxs]
              #powers2 = np.abs(sig2_fft)[pidxs]
              norm = len(sig)  # to normalize the power
              powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
              #powers2 = ((powers2/norm)**2)*(1./(sig2.std()**2))*2
              avg_array += powers
              #avg_array += powers2
          else:
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
          #plt.figure()
          #plt.loglog(freqs, powers2)
          #plt.ylim(10**-5, 10**0)    
          #plt.figure()
          #plt.loglog(freqs, powers)
          #plt.ylim(10**-5, 10**0)
        
        avg_array /= (2*n_segments)-1  # take the average of the segments
        """
        for w in range(1):
            sig = np.zeros((len(v_interp)/num_seg))
            #sig_test1 = sig[595:600]
            #print sig_test1
            percent_overlap = 99
            print len(sig)
            points_overlap = len(sig)*percent_overlap/100
            print points_overlap
            points_nonoverlap = (len(sig)-points_overlap)
            print points_nonoverlap
            for i in range(((len(v_interp)-len(sig))/points_nonoverlap)+1):  # maybe len(v_interp/#non-overlap points + 1)
                sig = v_interp[i*points_nonoverlap:(i*points_nonoverlap)+600] # i*#overlap points + freqs_size
                sig_fft = fftpack.fft(sig)
                powers = np.abs(sig_fft)[pidxs]
                norm = len(sig)  # to normalize the power
                powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
                avg_array += powers
            
                #plt.figure()
                #plt.loglog(freqs, powers)
                #plt.ylim(10**-5, 10**0)
                
            avg_array /= (len(v_interp)/points_nonoverlap)+1  # take the average of the segments
            #avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
            spectra_seg[ii-11][jj-11] = avg_array  # construct 3D array with averaged FFTs from each pixel
            
            fig = plt.figure(figsize=(15,12))
            plt.loglog(freqs, avg_array)
            #plt.loglog(freqs, avg_array)
            plt.ylim(10**-5, 10**0)
            plt.title('FFT 2-Hour Segment Averaging @ %i-percent overlap' % percent_overlap)
            
            
    # initialize arrays to hold temporary results for calculating arithmetic average (changed from geometric)
    temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
    p_avg = np.zeros((spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
    spectra_array = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
        
    
    ### calculate 3x3 pixel-box arithmetic average.  start at 1 and end 1 before to deal with edges.
    ## previously for geometric -- 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)
"""
for l in range(1,2):
#for l in range(1,25):
    #print l
    for m in range(1,2):
    #for m in range(1,25):

                   
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
        #spectra_array[l-1][m-1] = p_avg
        
        fig = plt.figure(figsize=(15,12))
        #plt.loglog(freqs, avg_array)
        plt.loglog(freqs, p_avg)
        plt.ylim(10**-5, 10**0)
        plt.title('FFT 2-Hour Segment Averaging @ %i-percent overlap' % percent_overlap)
        #plt.savefig('C:/Users/Brendan/Desktop/fft_2hr_segment_averaging_percent_overlap_%i.jpeg' % percent_overlap)
        #avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
                   
        #spectra_seg[ii][jj] = avg_array  # construct 3D array with averaged FFTs from each pixel
    
    
    # estimate time remaining and print to screen
"""    
"""    
# print estimated and total program time to screen 
print "Beginning Estimated time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec)
T_act = timer() - start
T_min3, T_sec3 = divmod(T_act, 60)
T_hr3, T_min3 = divmod(T_min3, 60)
print "Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3) 


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

return spectra_array
"""