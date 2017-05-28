# -*- coding: utf-8 -*-
"""
Created on Thu May 18 09:39:45 2017

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

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20130815'
wavelength = 171
num_seg = 1

#import accelerate  # switch on if computer has installed

from scipy import fftpack

"""
# in case need to add other descriptors to filename
flist_data = glob.glob('%s/DATA/Temp/%s/%i/*rebin1.npy' % (directory, date, wavelength)) 
flist_time = glob.glob('%s/DATA/Temp/%s/%i/*time.npy' % (directory, date, wavelength))
flist_exposure = glob.glob('%s/DATA/Temp/%s/%i/*exposure.npy' % (directory, date, wavelength))

part_name = '%s/DATA/Temp/%s/%i/' % (directory, date, wavelength)
len_data = len(flist_data[0])
len_time = len(flist_time[0])
len_exposure = len(flist_exposure[0])
len_part = len(part_name)

fname_data = data[0][len_part:len_data]
fname_time = data[0][len_part:len_time]
fname_exposure = data[0][len_part:len_exposure]

DATA = np.load('%s/DATA/Temp/%s/%i/%s' % (directory, date, wavelength, fname_data))

TIME = np.load('%s/DATA/Temp/%s/%i/%s' % (directory, date, wavelength, fname_time))

Ex = np.load('%s/DATA/Temp/%s/%i/%s' % (directory, date, wavelength, fname_exposure))
"""

DATA = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))

TIME = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))

Ex = np.load('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength))


print DATA.shape 

print "Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0])

    
## determine frequency values that FFT will evaluate
if wavelength == 1600 or wavelength == 1700:
  time_step = 24  # 24-second cadence for these wavelengths
else:
  time_step = 12  # 12-second cadence for the others
  
t_interp = np.linspace(0, TIME[len(TIME)-1], int((TIME[len(TIME)-1]/time_step)+1))  # interpolate onto default-cadence time-grid
  
n_segments = num_seg  # break data into 12 segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) / n_segments 

sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

reslt = (DATA.shape[0] == TIME.shape[0])
print "DATA and TIME array sizes match: %s" % reslt

pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_seg = np.zeros((9,len(freqs)))
spectra_std = np.zeros((len(freqs)))

print "length time-interp array = %i" % n
print "size for FFT to consider = %i" % freq_size
print "length of sample freq array = %i" % len(sample_freq)
print "length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs)


start = timer()
T1 = 0

for ii in range(163,166):
#for ii in range(0,5):

    for jj in range(245,248):
    #for jj in range(0,5):        
    
        #x1_box = 0+ii
        #x2_box = 2+ii  # if want to use median of more than 1x1 pixel box
        #y1_box = 0+jj
        #y2_box = 2+jj  # if want to use median of more than 1x1 pixel box
        
        """  # replace for-loop with direct assignment - 4x faster
        for k in range(0,DATA.shape[0]):
          #im=DATA[k]/(Ex[k])	  # get image + normalize by exposure time  (time went nuts?)
          im=DATA[k]
          #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])  # finds pixel-box median
          pixmed[k]= im[x1_box,y1_box]	# median  <-- use this
        """

        pixmed = DATA[:,ii,jj] / Ex  # extract timeseries + normalize by exposure time
        #pixmed = pixmed/Ex  # normalize by exposure time    
        

        
        v=pixmed  # 094/131/335 -- intensities too low to trim off negative values
        t=TIME
        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
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
        
                   
        spectra_seg[jj-245+(ii-163)*3] = powers  # construct 3D array with averaged FFTs from each pixel
        spectra_std = np.std(spectra_seg, axis=0)
        
        plt.figure(figsize=(15,15))
        plt.loglog(freqs,powers)
        plt.ylim(10**-6.5,10**0)
        
plt.figure(figsize=(15,15))
plt.loglog(freqs,spectra_std)

#np.save('C:/Users/Brendan/Desktop/spec_std.npy',spectra_std)
#np.save('C:/Users/Brendan/Desktop/spec_9_test.npy',spectra_seg)