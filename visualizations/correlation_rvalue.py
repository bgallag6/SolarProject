# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 21:18:22 2017

@author: Brendan
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)    

#m2 = np.load('F:/Users/Brendan/Desktop/SolarProject/Older Datafiles/20130626/20130626_1600_-450_-200i_-200_200j_m2.npy')
#spectra = np.load('F:/Users/Brendan/Desktop/SolarProject/Older Datafiles/20130626/20130626_193_-450_-200i_-200_200j_spectra.npy')
#param = np.load('F:/Users/Brendan/Desktop/SolarProject/Older Datafiles/20130626/20130626_193_-450_-200i_-200_200j_param.npy')

from scipy.stats.stats import pearsonr

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20130815'
wavelength = 193

spectra = np.load('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength))
param = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))

if wavelength == 1600:
    time_step = 24
else:
    time_step = 12  # add as argument, or leave in as constant?

if date == '20120923':   
    n_segments = 3  # break data into 12 segments of equal length
else:
    n_segments = 6
    
t_interp = np.linspace(0, 3600*2*n_segments, 3600*2*n_segments/time_step)

n = len(t_interp)
rem = n % n_segments  # n_segments is argument in module call (default to 1?)
freq_size = (n - rem) / n_segments

### determine frequency values that FFT will evaluate
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)    
freqs = sample_freq[pidxs]

m2 = np.zeros_like((spectra))

for i in range(spectra.shape[0]):
    for j in range(spectra.shape[1]):
        m2[i][j] = GaussPowerBase(freqs, param[0][i][j], param[1][i][j], param[2][i][j], param[3][i][j], param[4][i][j], param[5][i][j]) 
    

    


r = np.zeros((m2.shape[0], m2.shape[1]))

for i in range(m2.shape[0]):
    for j in range(m2.shape[1]):
        r_temp = pearsonr(m2[i][j], spectra[i][j])
        r[i][j] = r_temp[0]
        
    


if r.shape[1] > r.shape[0]:
    aspect_ratio = float(r.shape[1]) / float(r.shape[0])
    fig_height = 10
    fig_width = 10*aspect_ratio
    
else:
    aspect_ratio = float(r.shape[0]) / float(r.shape[1])
    #print aspect_ratio
    #fig_width = 10
    fig_width = 10+2  # works better for 20130626 (with no x/y labels)
    fig_height = 10*aspect_ratio  # works better for 20130626
#"""

#fig_width = 10+2  # works better for 20130626 (with no x/y labels)
#fig_width = 10+3  # works better for 20130626 (with x/y labels)
#fig_height = 10  # works better for 20130626

font_size = 23
    
#fig = plt.figure(figsize=(13,9))
fig = plt.figure(figsize=(fig_width,fig_height))
ax = plt.gca()  # get current axis -- to set colorbar 
#plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
plt.title('Correlation: r-Value', y = 1.02, fontsize=font_size)  # no date / wavelength

r = np.nan_to_num(r)
h_min = np.percentile(r,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
h_max = np.percentile(r,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
#cmap = 'jet'
cmap = cm.get_cmap('jet', 10)

# specify colorbar ticks to be at boundaries of segments
h_range = np.abs(h_max-h_min)
h_step = h_range / 10.
#h1 = h_min + (h_step/2.)
c_ticks = np.zeros((11))  # use 10 if centers of segments
for h in range(11):
    #c_ticks[h] = h1 + h_step*h  # use for centers of segments
    c_ticks[h] = h_min + h_step*h 
    
im = ax.imshow(np.flipud(r), cmap = cmap, vmin=h_min, vmax=h_max)
#im = ax.imshow(h_map[i], cmap = cmap, vmin=h_min, vmax=h_max)
#plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
#plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
#plt.xticks(x_ticks,fontsize=font_size)
#plt.yticks(y_ticks,fontsize=font_size)
#plt.xticks(x_ticks,x_ind,fontsize=font_size)
#plt.yticks(y_ticks,y_ind,fontsize=font_size)
ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax, format='%0.2f')
#cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
#cbar.set_ticks(np.round(c_ticks,8))  # 8 for slope (or might as well be for all, format separately)
cbar.set_ticks(c_ticks)  # 8 for slope (or might as well be for all, format separately)
#plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.jpeg' % (directory, date, wavelength, date, wavelength, names[i]))
#plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf')

