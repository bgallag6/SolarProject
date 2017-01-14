# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 11:29:02 2016

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

Y = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20130530_1600_0_296_0_634_numpy.npy')

R = np.zeros((Y.shape[1],Y.shape[2]))

for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            R[i][j] = Y[1][i][j]
            
for i in range(10,20):
    for j in range(10,20):
        R[i][j] = np.nan
            
fig, ax = plt.subplots()
vmin_initial = 0
vmax_initial = 3

ax1 = plt.subplot2grid((1,4),(0, 0), colspan=2, rowspan=1)
plt.subplots_adjust(bottom=0.25)
ax1.set_xlim(0, Y.shape[2]-1)
ax1.set_ylim(0, Y.shape[1]-1)  
im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # gauss amp

cbar = plt.colorbar(im,fraction=0.04, pad=0.01)
cbar.set_label("Gaussian Amplitude", size=20, labelpad=10)
cbar.set_label("Slope Value", size=20, labelpad=10)
cbar.set_label("$\chi_r^2$", size=20, labelpad=10)
cbar.ax.tick_params(labelsize=17, pad=1)    

ax2 = plt.subplot2grid((1,4),(0, 2), colspan=2, rowspan=1)

axcolor = 'white'
ax_vmin = plt.axes([0.05, 0.1, 0.45, 0.04], axisbg=axcolor)
s_vmin = Slider(ax_vmin, 'vmin', 0.1, 3., valinit=vmin_initial)




def update_min(val):
    vmin = s_vmin.val
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            R[i][j] = Y[1][i][j]
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            if R[i][j] < vmin:
                R[i][j] = np.nan
    #im = ax1.imshow(R, cmap='jet', interpolation='nearest', picker=True) 
    ax1 = plt.subplot2grid((1,4),(0, 0), colspan=2, rowspan=1)
    plt.subplots_adjust(bottom=0.25)
    ax1.set_xlim(0, Y.shape[2]-1)
    ax1.set_ylim(0, Y.shape[1]-1)  
    #im, = ([ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)])  # gauss amp
    #cbar, = ([plt.colorbar(im,fraction=0.04, pad=0.01)])
    im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # gauss amp
    cbar = plt.colorbar(im,fraction=0.04, pad=0.01)
    cbar.set_label("Gaussian Amplitude", size=20, labelpad=10)
    cbar.set_label("Slope Value", size=20, labelpad=10)
    cbar.set_label("$\chi_r^2$", size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=1)    
    
def update_max(val):
    vmax = s_vmax.val
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            R[i][j] = Y[1][i][j]
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            if R[i][j] > vmax:
                R[i][j] = np.nan
    #im = ax1.imshow(R, cmap='jet', interpolation='nearest', picker=True) 
    ax1 = plt.subplot2grid((1,4),(0, 0), colspan=2, rowspan=1)
    plt.subplots_adjust(bottom=0.25)
    ax1.set_xlim(0, Y.shape[2]-1)
    ax1.set_ylim(0, Y.shape[1]-1)  
    #im, = ([ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)])  # gauss amp
    #cbar, = ([plt.colorbar(im,fraction=0.04, pad=0.01)])
    im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # gauss amp
    cbar = plt.colorbar(im,fraction=0.04, pad=0.01)
    cbar.set_label("Gaussian Amplitude", size=20, labelpad=10)
    cbar.set_label("Slope Value", size=20, labelpad=10)
    cbar.set_label("$\chi_r^2$", size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=1)    
    


ax_vmax = plt.axes([0.05, 0.05, 0.45, 0.04], axisbg=axcolor)
s_vmax = Slider(ax_vmax, 'vmax', 0.1, 3., valinit=vmax_initial)
    #fig.canvas.draw_idle()
s_vmin.on_changed(update_min)
s_vmax.on_changed(update_max)


resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    s_vmin.reset()
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            R[i][j] = Y[1][i][j]
    ax1.imshow(R)    
button.on_clicked(reset)

plt.show()



"""
Y = np.load('F:/SDO/M2_Spectra_Params/param_20130530_1600_0_296_0_634_numpy.npy')


R = Y[1]

fig = plt.figure()            
plt.imshow(R)

v = 1.8

for i in range(R.shape[0]):
    for j in range(R.shape[1]):
        if R[i][j] < v:
            R[i][j] = np.NaN

fig = plt.figure()            
plt.imshow(R)
"""