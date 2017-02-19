# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 09:26:53 2016

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from matplotlib.widgets import  RectangleSelector
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
from scipy.interpolate import interp1d
from scipy import signal
import numpy as np
import scipy.misc
import astropy.units as u
import h5py
from scipy import fftpack

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)


    
DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_data_rebin1.npy')
TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_time.npy')
EXPOSURE = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_exposure.npy')

pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values


t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters

for ii in range(139,140):
    
    for jj in range(287,288):
        
        #print DATA[0][i][j]  # level / y / x
        
        x1_box = 0+ii
        x2_box = 1+ii
        y1_box = 0+jj
        y2_box = 1+jj
        
        for k in range(0,DATA.shape[0]):
                    im=DATA[k]												# get image
                    pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])	# median
        
        pixmed /= EXPOSURE           
                    
        # The derotation introduces some bad data towards the end of the sequence. This trims that off
        bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
        last_good_pos = bad - 1				# retain only data before the <=zero
        
        # Get t and v data
        v=pixmed[0:last_good_pos]		
        t=TIME[0:last_good_pos]
        
        
        
        #plt.plot(t,v)
    
        
        #t_interp = [l for l in range(0,TIME[len(TIME)-1])]
        v_interp = np.interp(t_interp,t,v)
        t_interp = np.array(t_interp)
        t_interp = t_interp.astype(float)
        
        
        data = v_interp
        avg_array = np.zeros((299))

        n_segments = 6  # break data into 12 segments of equal length
        n = len(data)
        r = n % n_segments
        data = data[0:len(data)-r]   # trim timeseries to be integer multiple of n_segments
        split = np.split(data, n_segments)
        for i in range(0,n_segments):

             ### FFT - Non-Detrended
             ### http://www.scipy-lectures.org/intro/scipy.html
             time_step = 12
             sig = split[i]
             sample_freq = fftpack.fftfreq(sig.size, d=time_step)
             sig_fft = fftpack.fft(sig)
             pidxs = np.where(sample_freq > 0)
             freqss = sample_freq[pidxs]
             powers = np.abs(sig_fft)[pidxs]
             norm = len(sig)  # to normalize the power
             powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
             avg_array += powers
            

             plt.figure(figsize=(20,15))
             plt.suptitle('Timeseries w/ FFT - 6 Segment Averaging', fontsize=20, fontweight='bold')
                            
             ax1 = plt.subplot2grid((16,3), (0,0), rowspan=4, colspan=3) 
             plt.plot(t_interp, v_interp, 'k')
             ax1.set_xlim([0, t_interp[len(t_interp)-1]])
             ax1.axvline(7196, color='r', linewidth=3)
             ax1.axvline(14392, color='r', linewidth=3)
             ax1.axvline(21588, color='r', linewidth=3)
             ax1.axvline(28784, color='r', linewidth=3)
             ax1.axvline(35980, color='r', linewidth=3)
             ax1.text(500 + (7196*i), 900, '%i' % (i+1), color='r', fontsize=20)
             plt.title('Timeseries', fontsize=20)
             plt.xlabel('Time [Seconds]', fontsize=15)
             plt.ylabel('Power', fontsize=15)
                                 
             ax2 = plt.subplot2grid((16,3), (6,0), rowspan=10, colspan=3)
             plt.loglog(freqss, powers, 'k')
             ax2.set_xlim([10**-4,10**-1])
             ax2.set_ylim([10**-6,10**0])
             plt.title('Fast Fourier Transform', fontsize=20)
             plt.xlabel('Frequency [Hz]', fontsize=15)
             plt.ylabel('Power', fontsize=15)
  
             #plt.savefig('C:/Users/Brendan/Desktop/PHYS 326/6_segment_averaging_%i.jpeg' % i)
             #plt.close()
  
          
        avg_array /= n_segments
        plt.figure(figsize=(20,15))
        plt.suptitle('Timeseries w/ FFT - 6 Segment Averaging', fontsize=20, fontweight='bold')
                   
        ax1 = plt.subplot2grid((16,3), (0,0), rowspan=4, colspan=3) 
        plt.plot(t_interp, v_interp, 'k')
        ax1.set_xlim([0, t_interp[len(t_interp)-1]])
        plt.title('Timeseries', fontsize=20)
        plt.xlabel('Time [Seconds]', fontsize=15)
        plt.ylabel('Power', fontsize=15)
                             
        ax2 = plt.subplot2grid((16,3), (6,0), rowspan=10, colspan=3)
        plt.loglog(freqss, avg_array, 'k')
        ax2.set_xlim([10**-4,10**-1])
        ax2.set_ylim([10**-6,10**0])
        plt.title('Fast Fourier Transform', fontsize=20)
        plt.xlabel('Frequency [Hz]', fontsize=15)
        plt.ylabel('Power', fontsize=15)
  
        #plt.savefig('C:/Users/Brendan/Desktop/PHYS 326/6_segment_averaging_full.jpeg')
        #plt.close()
