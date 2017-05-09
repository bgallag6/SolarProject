# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:47:41 2017

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
import accelerate  # switch on if computer has installed
from timeit import default_timer as timer

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)


    
#DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193/derotated.npy')
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193/time.npy')
#EXPOSURE = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193/exposure.npy')
DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130626/PCBs/171/derotated.npy')
TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130626/PCBs/171/time.npy')
EXPOSURE = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130626/PCBs/171/exposure.npy')

wavelength = 171




t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters

n_segments = 6  # break data into 12 segments of equal length  # change to [hours]
n = len(t_interp)
r = n % n_segments
freq_size = (n - r) / n_segments

if wavelength == 1600 or wavelength == 1700:
  time_step = 24  # 24 second cadence for these wavelengths
else:
  time_step = 12  # 12 second cadence for the others
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values

spec_array = np.zeros((11,DATA.shape[1],DATA.shape[2],len(freqs)))

start = timer()
T1 = 0

for ii in range(spec_array.shape[1]):
#for ii in range(139,140):
   
    for jj in range(spec_array.shape[2]):
    #for jj in range(287,300):
        
        #print DATA[0][i][j]  # level / y / x
        
        x1_box = 0+ii
        #x2_box = 1+ii
        y1_box = 0+jj
        #y2_box = 1+jj
        
        for k in range(0,DATA.shape[0]):
                    im=DATA[k]												# get image
                    #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])	# median
                    pixmed[k]=im[x1_box,y1_box]	# median
        
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
        avg_smooth = np.zeros((299))

        n_segments = 6  # break data into 12 segments of equal length  # change to [hours]
        n = len(data)
        r = n % n_segments
        data = data[0:len(data)-r]   # trim timeseries to be integer multiple of n_segments
        split = np.split(data, n_segments)
        time_step = 12  # add as argument in function call, or leave in as constant?
        freq_size = (n - r) / n_segments
        sample_freq = fftpack.fftfreq(freq_size, d=time_step)
        pidxs = np.where(sample_freq > 0)
        freqs = sample_freq[pidxs]

        
        #for w in range(1):
        sig = np.zeros((len(v_interp)/n_segments))
        #sig_test1 = sig[595:600]
        #print sig_test1
        percent_overlap = 50
        #print len(sig)
        points_overlap = len(sig)*percent_overlap/100
        #print points_overlap
        points_nonoverlap = (len(sig)-points_overlap)
        #print points_nonoverlap
        for i in range(((len(v_interp)-len(sig))/points_nonoverlap)+1):  # maybe len(v_interp/#non-overlap points + 1)
            sig = v_interp[i*points_nonoverlap:(i*points_nonoverlap)+600] # i*#overlap points + freqs_size
            #print (i*points_nonoverlap), (i*points_nonoverlap)+600
            #sig_fft = fftpack.fft(sig)
            sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
            powers = np.abs(sig_fft)[pidxs]
            norm = len(sig)  # to normalize the power
            powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
            avg_array += powers
                  
            orig_powerspec = powers
            # Now Smooth the Power Spectra
            orig = orig_powerspec.copy()		# make a copy            
            window_len=7   # Keep this an odd number            
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
            avg_smooth += y2
            
            avg_raw = np.copy(avg_array)  # take the average of the segments
            avg_raw = avg_raw / (i+1)
            
            avg_smooth2 = np.copy(avg_smooth)
            avg_smooth2 = avg_smooth2 / (i+1)
            
            #spec_array[0][i] = powers
            spec_array[i][ii][jj] = y2
            
            """
            plt.figure(figsize=(20,15))
            plt.suptitle('Timeseries w/ FFT - 6 Segment Averaging', fontsize=20, fontweight='bold')
                            
            ax1 = plt.subplot2grid((16,4), (0,0), rowspan=4, colspan=4) 
            plt.plot(t_interp, v_interp, 'k')
            ax1.set_xlim([0, t_interp[len(t_interp)-1]])
            ax1.axvline(t_interp[points_nonoverlap*i], color='r', linewidth=3)
            ax1.axvline(t_interp[points_nonoverlap*i+600-1], color='r', linewidth=3)
            #ax1.axvline(21588, color='r', linewidth=3)
            #ax1.axvline(28784, color='r', linewidth=3)
            #ax1.axvline(35980, color='r', linewidth=3)
            #ax1.text(500 + (7196*i), 900, '%i' % (i+1), color='r', fontsize=20)
            plt.title('Timeseries', fontsize=20)
            plt.xlabel('Time [Seconds]', fontsize=15)
            plt.ylabel('Power', fontsize=15)
            
            ax2 = plt.subplot2grid((16,4), (6,0), rowspan=10, colspan=2)
            #plt.loglog(freqs[2:], powers[2:], 'k')
            #plt.loglog(freqs[2:], y2[1:len(y2)-1], 'r')
            plt.loglog(freqs, powers, 'k')
            plt.loglog(freqs[5:], y2[3:len(y2)-2], 'r')
            ax2.set_xlim([10**-4,10**-1])
            ax2.set_ylim([10**-6,10**0])
            plt.title('Raw Spectra', fontsize=20)
            plt.xlabel('Frequency [Hz]', fontsize=15)
            plt.ylabel('Power', fontsize=15)
            
            ax3 = plt.subplot2grid((16,4), (6,2), rowspan=10, colspan=2)
            #plt.loglog(freqs[2:], avg_raw[2:], 'k')
            #plt.loglog(freqs[2:], avg_smooth2[1:len(avg_smooth2)-1], 'r')
            plt.loglog(freqs, avg_raw, 'k')
            plt.loglog(freqs[5:], avg_smooth2[3:len(avg_smooth2)-2], 'r')
            ax3.set_xlim([10**-4,10**-1])
            ax3.set_ylim([10**-6,10**0])
            plt.title('Running Averages', fontsize=20)
            plt.xlabel('Frequency [Hz]', fontsize=15)
            plt.ylabel('Power', fontsize=15)
            
            #plt.savefig('C:/Users/Brendan/Desktop/fft_overlap2/segment_%i.jpeg' % i)
            #plt.figure()
            #plt.loglog(freqs, powers)
            #plt.ylim(10**-5, 10**0)
            """
        
        #np.save('C:/Users/Brendan/Desktop/fft_overlap2/array_75pct.npy', spec_array)
        #avg_array /= (len(v_interp)/points_nonoverlap)+1  # take the average of the segments
        #avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
        #spectra_seg[ii-11][jj-11] = avg_array  # construct 3D array with averaged FFTs from each pixel        
        # estimate time remaining and print to screen
    T = timer()
    T2 = T - T1
    if ii == 0:
        T_init = T - start
        T_est = T_init*(spec_array.shape[1])  
        T_min, T_sec = divmod(T_est, 60)
        T_hr, T_min = divmod(T_min, 60)
        #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est)
        print "Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spec_array.shape[1], T_hr, T_min, T_sec)
    else:
        T_est2 = T2*(spec_array.shape[1]-ii)
        T_min2, T_sec2 = divmod(T_est2, 60)
        T_hr2, T_min2 = divmod(T_min2, 60)
        #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est2)
        print "Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spec_array.shape[1], T_hr2, T_min2, T_sec2)
    T1 = T
        
"""
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
        """
"""  
        #avg_array /= n_segments
        plt.figure(figsize=(20,15))
        plt.suptitle('Timeseries w/ FFT - 6 Segment Averaging', fontsize=20, fontweight='bold')
                   
        ax1 = plt.subplot2grid((16,3), (0,0), rowspan=4, colspan=3) 
        plt.plot(t_interp, v_interp, 'k')
        ax1.set_xlim([0, t_interp[len(t_interp)-1]])
        plt.title('Timeseries', fontsize=20)
        plt.xlabel('Time [Seconds]', fontsize=15)
        plt.ylabel('Power', fontsize=15)
                             
        ax2 = plt.subplot2grid((16,3), (6,0), rowspan=10, colspan=3)
        #plt.loglog(freqss, avg_array, 'k')
        plt.loglog(freqs, avg_array, 'k')
        ax2.set_xlim([10**-4,10**-1])
        ax2.set_ylim([10**-6,10**0])
        plt.title('Fast Fourier Transform', fontsize=20)
        plt.xlabel('Frequency [Hz]', fontsize=15)
        plt.ylabel('Power', fontsize=15)
        
        #plt.savefig('C:/Users/Brendan/Desktop/fft_overlap2/segment_averaged.jpeg')
        #plt.close()
        """

np.save('C:/Users/Brendan/Desktop/fft_overlap2/20130626_171_PCB_50pct_2hr.npy',spec_array)