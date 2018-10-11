# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 11:29:02 2016

@author: Brendan
"""

# update 1/13: took out spectra plot on side - not sure if want

# improved formatting a bit - copied over heatmap tool colorbar setup + percentiles

# possibly add a text input for vmin and vmax?

# add the different parameter heatmaps as buttons at the top 
# - will require a bunch of generalization - for the vmin/vmax 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from pylab import *
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import figure, show

directory = 'F:'
date = '20130626'
wavelength = 304


h_map0 = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
#h_map = np.copy(h_map0)
#h_map = np.load('C:/Users/Brendan/Desktop/SDO/param_20130530_193_2300_2600_2200_3000.npy')

#wavelength = 211
#year = 2012
#month = 9
#day = 23


def update_min(val):
    vmin = s_vmin.val
#    print vmin
    h_map = np.copy(h_map0)
    #for i in range(R.shape[0]):
    #    for j in range(R.shape[1]):
    #        R[i][j] = h_map[1][i][j]
    #for i in range(R.shape[0]):
    #    for j in range(R.shape[1]):
    #        if R[i][j] < vmin:
    #            R[i][j] = np.nan
    R = h_map[p]
    R[R<vmin] = np.NaN
 
    ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
    plt.subplots_adjust(bottom=0.25)
    plt.subplots_adjust(left=0.05)
    ax1.set_xlim(0, h_map.shape[2]-1)
    ax1.set_ylim(0, h_map.shape[1]-1)   
    ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
    #im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # before revision
    im = ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
    cbar.ax.tick_params(labelsize=13, pad=3)
    
def update_max(val):
    vmax = s_vmax.val
#    print vmax
    h_map = np.copy(h_map0)
    #for i in range(R.shape[0]):
    #    for j in range(R.shape[1]):
    #        R[i][j] = h_map[1][i][j]
    #for i in range(R.shape[0]):
    #    for j in range(R.shape[1]):
    #        if R[i][j] > vmax:
    #            R[i][j] = np.nan
    R = h_map[p]
    R[R>vmax] = np.NaN

    ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
    plt.subplots_adjust(bottom=0.25)
    plt.subplots_adjust(left=0.05)
    ax1.set_xlim(0, h_map.shape[2]-1)
    ax1.set_ylim(0, h_map.shape[1]-1)  
    ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
    #im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # before revision
    im = ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
    cbar.ax.tick_params(labelsize=13, pad=3)    
    

def reset(event):
    s_vmin.reset()
    s_vmax.reset()
    h_map = np.copy(h_map0)
    R = h_map[p]
    ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)    
    
    

p = 1

#R = np.zeros((h_map.shape[1],h_map.shape[2]))

#date_title = '%i-%02i-%02i' % (year,month,day)
date_title = date
        
# create list of titles and colorbar names for display on the figures
titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width', '$/chi^2$', 'Visual Image - Averaged']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '$/chi^2$', 'Intensity']

h_map = np.copy(h_map0)
R = h_map[p]
#for i in range(R.shape[0]):
#        for j in range(R.shape[1]):
#            R[i][j] = h_map[1][i][j]
            

# create figure with heatmap and spectra side-by-side subplots
fig1 = plt.figure(figsize=(20,10))
ax1 = plt.gca()
ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
plt.subplots_adjust(bottom=0.25)
plt.subplots_adjust(left=0.05)
ax1.set_xlim(0, h_map.shape[2]-1)
ax1.set_ylim(0, h_map.shape[1]-1)  
ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[p]), y = 1.01, fontsize=17)
                  


param = h_map[p]  # set initial heatmap to power law index     
h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
im = ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)

#fig, ax = plt.subplots()
vmin_initial = h_min
vmax_initial = h_max
#im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # gauss amp

# design colorbar for heatmaps
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
cbar.ax.tick_params(labelsize=13, pad=3)   
        

#ax2 = plt.subplot2grid((1,11),(0, 6), colspan=5, rowspan=1)
#ax2.loglog()
#ax2.set_xlim(10**-4.5, 10**-1.3)
#ax2.set_ylim(10**-5, 10**0)  

sl_min = [10**-9, 0.3, 10**-7, 10**-5, -6.4, 0.05]
sl_max = [10**-4, 4.0, 10**-2, 10**-1, -4.6, 0.8]

axcolor = 'white'
ax_vmin = plt.axes([0.2, 0.13, 0.45, 0.04], axisbg=axcolor)
s_vmin = Slider(ax_vmin, 'vmin', sl_min[p], sl_max[p], valinit=vmin_initial)

ax_vmax = plt.axes([0.2, 0.08, 0.45, 0.04], axisbg=axcolor)
s_vmax = Slider(ax_vmax, 'vmax', sl_min[p], sl_max[p], valinit=vmax_initial)
    #fig.canvas.draw_idle()
s_vmin.on_changed(update_min)
s_vmax.on_changed(update_max)


resetax = plt.axes([0.72, 0.105, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

button.on_clicked(reset)

plt.show()