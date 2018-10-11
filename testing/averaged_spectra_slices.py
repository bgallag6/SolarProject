# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 07:33:30 2018

@author: Brendan
"""

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import h5py
from scipy import fftpack  # doesnt work in module when called here???
import matplotlib.pylab as plt
from astropy.convolution import convolve, Box1DKernel
from matplotlib import cm
from numpy.random import randn
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.ticker import LogFormatterMathtext
from timeit import default_timer as timer
from scipy.stats import f
import matplotlib.patches as patches
from scipy.stats import f as ff
from scipy.stats.stats import pearsonr
import matplotlib   

   
directory = 'F:'
date = 20130626
wavelength = 171
  
cube = np.load('%s/DATA/Temp/%s/PCBs/%i/derotated.npy' % (directory, date, wavelength))

TIME = np.load('%s/DATA/Temp/%s/PCBs/%i/time.npy' % (directory, date, wavelength))
exposure = np.load('%s/DATA/Temp/%s/PCBs/%i/exposure.npy' % (directory, date, wavelength))

# determine frequency values that FFT will evaluate
if wavelength in [1600,1700]:
    time_step = 24  # add as argument in function call, or leave in as constant?
else:
    time_step = 12

t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/time_step)+1)  # interpolate onto default-cadence time-grid
    
n_segments = 6  # break data into # segments of equal length
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) // n_segments

sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

DATA = cube
    
    
pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))

ii = 235
jj = 238

pixmed = DATA[:,ii,jj] / exposure  # extract timeseries + normalize by exposure time   

v_interp = np.interp(t_interp,TIME,pixmed)  # interpolate pixel-intensity values onto specified time grid

data = v_interp

avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers

data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
split = np.split(data, n_segments)  # create split array for each segment
 
plt.rcParams["font.family"] = "Times New Roman"
font_size = 27

  
for i in range(0,n_segments):               
    
  ## perform Fast Fourier Transform on each segment       
  sig = split[i]
  sig_fft = fftpack.fft(sig)
  powers = np.abs(sig_fft)[pidxs]
  norm = len(sig)
  powers = ((powers/norm)**2)*(1./(sig.std()**2))*2   # normalize the power
  avg_array += powers  
  
  fig = plt.figure(figsize=(12,10))
  ax = plt.gca()  # get current axis -- to set colorbar 
  plt.title(r'2-hour slices: %i of %i' % ((i+1), 6), y = 1.01, fontsize=font_size)
  plt.ylim((10**-4.7,10**0))
  plt.xlim((10**-4.,10**-1.3))
  plt.xticks(fontsize=font_size)
  plt.yticks(fontsize=font_size)
  ax.tick_params(axis='both', which='major', pad=10)
  plt.loglog(freqs,powers,'k', label='Averaged Spectrum')

  plt.xlabel(r'Frequency [Hz]', fontsize=font_size, labelpad=10)
  plt.ylabel(r'Power', fontsize=font_size, labelpad=10)
  plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
  plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
  #plt.savefig('C:/Users/Brendan/Desktop/%i_of_%i.pdf' % ((i+1), 6), format='pdf', bbox_inches='tight')

avg_array /= n_segments  # take the average of the segments

avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
  
     
fig = plt.figure(figsize=(12,10))
ax = plt.gca()  # get current axis -- to set colorbar 
plt.title(r'Averaged Spectrum', y = 1.01, fontsize=font_size)
plt.ylim((10**-4.7,10**0))
plt.xlim((10**-4.,10**-1.3))
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax.tick_params(axis='both', which='major', pad=10)
plt.loglog(freqs,avg_array,'k', label='Averaged Spectrum')


plt.xlabel(r'Frequency [Hz]', fontsize=font_size, labelpad=10)
plt.ylabel(r'Power', fontsize=font_size, labelpad=10)
plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
plt.savefig('C:/Users/Brendan/Desktop/full.pdf', format='pdf', bbox_inches='tight')