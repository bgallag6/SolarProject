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
from mpl_toolkits.axes_grid1 import make_axes_locatable

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)


    
DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193_small/derotated.npy')
TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193_small/time.npy')
EXPOSURE = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130530/193_small/exposure.npy')

VISUAL = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130530/193/193_original_region/visual.npy')
V_AVG = VISUAL[0,132:145,279:292]

v_min = np.percentile(V_AVG,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
v_max = np.percentile(V_AVG,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%) 

plt.figure(figsize=(12,9))

ax = plt.gca()
#plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
plt.title(r'193$\AA$ - Visual: Average', y = 1.02, fontsize=20)  # no date / wavelength
#plt.title('3x3 Pixel-Box Averaging', fontsize=30, fontweight='bold', y=1.02)
im = ax.imshow(V_AVG, cmap='sdoaia193', vmin = v_min, vmax = v_max)
plt.xlabel('X-Position [Pixels]', fontsize=15, labelpad=8)
plt.ylabel('Y-Position [Pixels]', fontsize=15, labelpad=8)
ax.tick_params(axis='both', which='major', pad=8)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
#cbar.set_label('Intensity', size=20, labelpad=10)
cbar.ax.tick_params(labelsize=15, pad=5) 
rect = patches.Rectangle((6,5), 1, 1, color='red', alpha=.80, fill=True)
ax.add_patch(rect)
rect = patches.Rectangle((7,6), 1, 1, color='red', fill=False)
ax.add_patch(rect)  
rect = patches.Rectangle((5,4), 1, 1, color='red', fill=False)
ax.add_patch(rect)
rect = patches.Rectangle((5,6), 1, 1, color='red', fill=False)
ax.add_patch(rect)
rect = patches.Rectangle((7,5), 1, 1, color='red', fill=False)
ax.add_patch(rect)  
rect = patches.Rectangle((7,4), 1, 1, color='red', fill=False)
ax.add_patch(rect)
rect = patches.Rectangle((5,5), 1, 1, color='red', fill=False)
ax.add_patch(rect)
rect = patches.Rectangle((6,4), 1, 1, color='red', fill=False)
ax.add_patch(rect)  
rect = patches.Rectangle((6,6), 1, 1, color='red', fill=False)
ax.add_patch(rect)
plt.savefig('C:/Users/Brendan/Desktop/3x3_visualB.jpeg')



pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
spectra_seg = np.zeros((3,3,299))

t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters


"""
for ii in range(139,140):
    for jj in range(288,289):

#for ii in range(138,141):
#    for jj in range(287,290):
        
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
             ax1.text(500 + (7196*i), 800, '%i' % (i+1), color='r', fontsize=20)
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
        plt.suptitle('6-Segment Averaging', fontsize=23, fontweight='bold')
                   
        ax1 = plt.subplot2grid((16,3), (0,0), rowspan=4, colspan=3) 
        plt.plot(t_interp, v_interp, 'k')
        ax1.set_xlim([0, t_interp[len(t_interp)-1]])
        ax1.axvline(7196, color='r', linewidth=3)
        ax1.axvline(14392, color='r', linewidth=3)
        ax1.axvline(21588, color='r', linewidth=3)
        ax1.axvline(28784, color='r', linewidth=3)
        ax1.axvline(35980, color='r', linewidth=3)
        ax1.text(500 + (7196*0), 900, '%i' % (0+1), color='r', fontsize=20)
        ax1.text(500 + (7196*1), 900, '%i' % (1+1), color='r', fontsize=20)
        ax1.text(500 + (7196*2), 900, '%i' % (2+1), color='r', fontsize=20)
        ax1.text(500 + (7196*3), 900, '%i' % (3+1), color='r', fontsize=20)
        ax1.text(500 + (7196*4), 900, '%i' % (4+1), color='r', fontsize=20)
        ax1.text(500 + (7196*5), 900, '%i' % (5+1), color='r', fontsize=20)
        plt.title('Time Series', fontsize=20, y=1.01)
        plt.xlabel('Time [Seconds]', fontsize=20)
        plt.ylabel('Intensity', fontsize=20)
                             
        ax2 = plt.subplot2grid((16,3), (6,0), rowspan=10, colspan=3)
        plt.loglog(freqss, powers, 'r', label = 'Single Segment')
        plt.loglog(freqss, avg_array, 'k', label = '6-Segment Averaged')
        ax2.set_xlim([10**-4,10**-1])
        ax2.set_ylim([10**-6,10**0])
        plt.title('FFT-Derived Power Spectra', fontsize=20, y=1.01)
        plt.xlabel('Frequency [Hz]', fontsize=20)
        plt.ylabel('Power', fontsize=20)
        legend = ax2.legend(loc='upper right', prop={'size':25}, labelspacing=0.35)
        for label in legend.get_lines():
            label.set_linewidth(2.5)  # the legend line width
  
        #plt.savefig('C:/Users/Brendan/Desktop/6_segment_averaging_full.jpeg')
        #plt.close()
        """  
"""  
        avg_array /= n_segments  # take the average of the segments
            
        avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array  
  
        spectra_seg[ii-138][jj-287] = avg_array  # construct 3D array with averaged FFTs from each pixel
        
        
# initialize arrays to hold temporary results for calculating arithmetic average (changed from geometric)
temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
p_avg = np.zeros((spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
spectra_array = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
    

### calculate 3x3 pixel-box arithmetic average.  start at 1 and end 1 before to deal with edges.
## previously for geometric -- 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

for l in range(1,spectra_seg.shape[0]-1):
#for l in range(1,25):
    #print l
    for m in range(1,spectra_seg.shape[1]-1):
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
        spectra_array[l-1][m-1] = p_avg
        
        
        plt.figure(figsize=(20,15))        
        plt.loglog(freqss, spectra_seg[l][m], 'r', label='Single Pixel')
        plt.loglog(freqss, p_avg, 'k', label='3x3 Averaged')
        plt.xlim([10**-4,10**-1])
        plt.ylim([10**-5,10**0])
        plt.title('3x3 Pixel-Box Averaging', fontsize=30, fontweight='bold', y=1.02)
        plt.xlabel('Frequency [Hz]', fontsize=20)
        plt.ylabel('Power', fontsize=20)
        legend = plt.legend(loc='upper right', prop={'size':30}, labelspacing=0.35)
        for label in legend.get_lines():
            label.set_linewidth(2.5)  # the legend line width
        plt.savefig('C:/Users/Brendan/Desktop/3x3_averaging.jpeg')
"""