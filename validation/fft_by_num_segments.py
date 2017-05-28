# -*- coding: utf-8 -*-
"""
Created on Thu May 18 11:18:22 2017

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
import matplotlib.pyplot as plt
from scipy import fftpack


directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20120909'  # was 20130815, 171 for other figures
wavelength = 1600

DATA = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))

TIME = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))

#Ex = np.load('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength))

#spec_std = np.load('C:/Users/Brendan/Desktop/spec_std.npy')
#spec_segs = np.load('C:/Users/Brendan/Desktop/spec_9_test.npy')

## determine frequency values that FFT will evaluate
if wavelength == 1600 or wavelength == 1700:
  #time_step = 24  # 24-second cadence for these wavelengths
  time_step = 12  # Jacks
else:
  time_step = 12  # 12-second cadence for the others
  
t_interp = np.linspace(0, TIME[len(TIME)-1], int((TIME[len(TIME)-1]/time_step)+1))  # interpolate onto default-cadence time-grid

#pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
pixmed=np.empty(DATA.shape[2])  # Initialize array to hold median pixel values


#pixmed = DATA[:,164,246] / Ex  # extract timeseries + normalize by exposure time
pixmed = DATA[10,10,:]  # Jacks

v=pixmed  # 094/131/335 -- intensities too low to trim off negative values
t=TIME
v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid


num_seg = np.array([1,3,6,12])   
    

for n in range(len(num_seg)):
 
    n_segments = num_seg[n]  # break data into 12 segments of equal length
    r = len(t_interp)
    rem = r % n_segments
    freq_size = (r - rem) / n_segments 
    
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    if n_segments == 1:
        freqs1 = freqs
    elif n_segments == 3:
        freqs3 = freqs  
    elif n_segments == 6:
        freqs6 = freqs
    elif n_segments == 12:
        freqs12 = freqs
           
    data = v_interp 
    data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
    split = np.split(data, n_segments)  # create split array for each segment
    avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
    avg_array2 = []
    #temp = np.zeros((len(freqs)))

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
      #temp += np.log10(powers)
    
    avg_array /= n_segments  # take the average of the segments 
    #temp /= n_segments
    #spec_geo = np.power(10,temp)
    
    if n_segments == 1:
        avg_array1 = avg_array
        #avg_array1 = spec_geo
    elif n_segments == 3:
        avg_array3 = avg_array/3.
        #avg_array3 = spec_geo/3.
    elif n_segments == 6:
        avg_array6 = avg_array/6.
        #avg_array6 = spec_geo/6.
    elif n_segments == 12:
        avg_array12 = avg_array/12.
        #avg_array12 = spec_geo/12.

#spectra_seg[jj-245+(ii-163)*3] = powers  # construct 3D array with averaged FFTs from each pixel
#spectra_std = np.std(spectra_seg, axis=0)

#plt.figure(figsize=(15,15))
#plt.loglog(freqs,avg_array)
#plt.ylim(10**-6.5,10**0)
        
plt.rcParams["font.family"] = "Times New Roman"
font_size = 27
    
plt.figure(figsize=(15,15))
ax = plt.gca()  
plt.title('Comparison of Temporal Averaging Methods \n Normalized by # of Segments | Jacks Data', y=1.01, fontsize=25)
plt.loglog(freqs1,avg_array1, 'k', linewidth=1.7, label='(1) 12-Hour Segment')
plt.loglog(freqs3,avg_array3, 'b', linewidth=1.7, label='(3) 4-Hour Segments')
plt.loglog(freqs6,avg_array6, 'g', linewidth=1.7, label='(6) 2-Hour Segments')
plt.loglog(freqs12,avg_array12, 'r', linewidth=1.7, label='(12) 1-Hour Segments')
plt.ylim(10**-6.5,10**0)
plt.xlim(10**-5.,10**-1.3)
plt.xticks(fontsize=font_size, fontname="Times New Roman")
plt.yticks(fontsize=font_size, fontname="Times New Roman")
plt.tick_params(axis='both', which='major', pad=10)
legend = plt.legend(loc='lower left', prop={'size':font_size}, labelspacing=0.35)
for label in legend.get_lines():
    label.set_linewidth(3.0)  # the legend line width
#plt.savefig('C:/Users/Brendan/Desktop/Averaging_Normalized_Jacks_300.pdf', format='pdf')

#np.save('C:/Users/Brendan/Desktop/spec_array1_jack.npy', avg_array1)
#np.save('C:/Users/Brendan/Desktop/spec_array3_jack.npy', avg_array3)
#np.save('C:/Users/Brendan/Desktop/spec_array6_jack.npy', avg_array6)
#np.save('C:/Users/Brendan/Desktop/spec_array12_jack.npy', avg_array12)