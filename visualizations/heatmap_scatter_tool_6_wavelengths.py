# -*- coding: utf-8 -*-
"""
Created on Tue May 23 17:14:18 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, ion, draw
import numpy as np
import random
from matplotlib.widgets import RectangleSelector
import matplotlib.patches as patches
from PIL import Image
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from matplotlib.widgets import Cursor
from pylab import *
from pylab import axvline
import sunpy
from sunpy.map import Map
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
from scipy import fftpack
from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib.pylab as plt
from astropy.convolution import convolve, Box1DKernel
from matplotlib import cm
from numpy.random import randn
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.ticker import LogFormatterMathtext
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib.pyplot import figure, show
from matplotlib.image import AxesImage
from matplotlib.widgets import Button
from matplotlib.widgets import RadioButtons
from matplotlib.pyplot import figure, show
import numpy as npy
from numpy.random import rand

    
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
    
class Index(object):
    ind = 0
    
       
    def coeff(self, event):
        param = h_map[0]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[0]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[0], size=15, labelpad=10)
        plt.draw()

    def index(self, event):
        param = h_map[1]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
        plt.draw()
        
    def roll(self, event):
        paramA = h_map[0]
        paramn = h_map[1]
        paramC = h_map[2]
        param = (paramC/paramA)**(1./paramn)
        param = np.nan_to_num(param)/60.
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[2]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[2], size=15, labelpad=10)
        plt.draw()
        
    def gauss_amp(self, event):
        param = h_map[3]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[3]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[3], size=15, labelpad=10)
        plt.draw()
        
    def gauss_loc(self, event):
        param = (1./(np.exp(h_map[4]))/60.)
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet_r', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[4]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[4], size=15, labelpad=10)
        plt.draw()
        
    def gauss_wid(self, event):
        param = h_map[5]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[5]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[5], size=15, labelpad=10)
        plt.draw()
        
    def fstat(self, event):
        param = h_map[8]
        NaN_replace = np.nan_to_num(param)  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='none', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[6]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[6], size=15, labelpad=10)
        plt.draw()
        
    def visual(self, event):
        param = vis[0]
        im = ax1.imshow(param, cmap='sdoaia%i' % wavelength, interpolation='nearest', picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[7]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[7], size=15, labelpad=10)
        plt.draw()
        
       
    def f1700(self, event):
        global t1700
        if t1700 == 0:
            t1700 = 1
        elif t1700 == 1:
            t1700 = 0
            plt.draw()                
        return t1700
        
    def f1600(self, event):
        global t1600
        if t1600 == 0:
            t1600 = 1
        elif t1600 == 1:
            t1600 = 0
            plt.draw()                
        return t1600
        
    def f304(self, event):
        global t304
        if t304 == 0:
            t304 = 1
        elif t304 == 1:
            t304 = 0
            plt.draw()                
        return t304
        
    def f171(self, event):
        global t171
        if t171 == 0:
            t171 = 1
        elif t171 == 1:
            t171 = 0
            plt.draw()                
        return t171
        
    def f193(self, event):
        global t193
        if t193 == 0:
            t193 = 1
        elif t193 == 1:
            t193 = 0
            plt.draw()                
        return t193
    
    def f211(self, event):
        global t211
        if t211 == 0:
            t211 = 1
        elif t211 == 1:
            t211 = 0
            plt.draw()                
        return t211

        
  

# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy, c, t1700, t1600, t304, t171, t193, t211
    ixx, iyy = event.xdata, event.ydata
    ax2.clear()
    plt.draw()
    print ('x = %d, y = %d' % ( ixx, iyy))  # print location of pixel
    ix = int(ixx)
    iy = int(iyy)
    
    #s1700 = np.zeros((spectra1700.shape[2]))
    s1600 = np.zeros((spectra1600.shape[2]))
    s304 = np.zeros((spectra304.shape[2]))
    s171 = np.zeros((spectra171.shape[2]))
    s193 = np.zeros((spectra193.shape[2]))
    s211 = np.zeros((spectra211.shape[2]))
    #m = np.zeros((spectra.shape[2]))
    #g = np.zeros((spectra.shape[2]))
    #pl = np.zeros((spectra.shape[2]))
    
    #s1700[:] = spectra1700[iy][ix][:]
    s1600[:] = spectra1600[iy][ix][:]
    s304[:] = spectra304[iy][ix][:]
    s171[:] = spectra171[iy][ix][:]
    s193[:] = spectra193[iy][ix][:]
    s211[:] = spectra211[iy][ix][:]
    #for i in range(0,spectra.shape[2]):
    #    s[i] = cpy_arr_spec[iy][ix][i]
    
    #m = GaussPowerBase(f_fit, param1[0][iy][ix], param1[1][iy][ix], param1[2][iy][ix], param1[3][iy][ix], param1[4][iy][ix], param1[5][iy][ix]) 
    #g = Gauss(f_fit, param1[3][iy][ix], param1[4][iy][ix], param1[5][iy][ix])
    #pl = PowerLaw(f_fit, param1[0][iy][ix], param1[1][iy][ix], param1[2][iy][ix])
    if t1700 == 1:
        ax2.loglog(f24, s1700, 'purple', label=r'1700$\AA$')
    if t1600 == 1:
        ax2.loglog(f24, s1600, 'green', label=r'1600$\AA$')
    if t304 == 1:
        ax2.loglog(f12, s304, 'red', label=r'304$\AA$')
    if t171 == 1:
        ax2.loglog(f12, s171, 'blue', label=r'171$\AA$')
    if t193 == 1:
        ax2.loglog(f12, s193, 'orange', label=r'193$\AA$')
    if t211 == 1:
        ax2.loglog(f12, s211, 'black', label=r'211$\AA$')
    
    #if toggle == 1: 
    #    ax2.loglog(f_fit1, s1, 'purple')
    #    ax2.loglog(f_fit1, s2, 'green')
    #ax2.loglog(f_fit, m, 'purple', label='M2 Combined')
    #ax2.loglog(f_fit, g, 'g--', label='Gaussian')
    #ax2.loglog(f_fit, pl, 'g', label='Power Law')
     
    ax2.set_xlim(10**-4.5, 10**-1.3)
    ax2.set_ylim(10**-5, 10**0)  
    axvline(x=0.00333,color='k',ls='dashed', label='5 minutes')
    axvline(x=0.00555,color='k',ls='dotted', label='3 minutes')
    ax2.set_title('Spectra Fit : Pixel (%ii , %ij)' % (iy, ix), fontsize=15)
    #plt.text(0.006, 10**-0.31, r'$A$ = {0:0.2e}'.format(h_map[0][iy][ix]), fontsize=23)
    #plt.text(0.0061, 10**-0.51, r'$n$ = {0:0.2f}'.format(h_map[1][iy][ix]), fontsize=23)
    #plt.text(0.006, 10**-0.73, r'$R$ = %0.1f [min]' % ((1./(h_map[2][iy][ix] / h_map[0][iy][ix])**(-1./ h_map[1][iy][ix]))/60.), fontsize=23)
    #plt.text(0.0061, 10**-0.95, r'$\alpha$ = {0:0.2e}'.format(h_map[3][iy][ix]), fontsize=23)
    #plt.text(0.0061, 10**-1.15, r'$\beta$ = {0:0.1f} [min]'.format((1./np.exp(h_map[4][iy][ix]))/60.), fontsize=23)
    #plt.text(0.0061, 10**-1.35, r'$\sigma$ = {0:0.3f}'.format(h_map[5][iy][ix]), fontsize=23)
    #plt.legend(loc='lower left', prop={'size':20})
    legend = ax2.legend(loc='lower left', prop={'size':25}, labelspacing=0.35)
    for label in legend.get_lines():
            label.set_linewidth(3.0)  # the legend line width
    plt.draw()

    return ix, iy
    
# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)    

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
        
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2) 

directory = 'D:/Users/Brendan/Desktop/SolarProject'
date = '20130626'
wavelength1700 = 1700
wavelength1600 = 1600
wavelength304 = 304
wavelength171 = 171
wavelength193 = 193
wavelength211 = 211

global spectra1700
global spectra1600
global spectra304
global spectra171
global spectra193
global spectra211

#cube_shape1700 = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength1700))
#spectra1700 = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength1700), dtype='float64', mode='r', shape=(cube_shape1700[0], cube_shape1700[1], cube_shape1700[2]))
cube_shape1600 = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength1600))
spectra1600 = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength1600), dtype='float64', mode='r', shape=(cube_shape1600[0], cube_shape1600[1], cube_shape1600[2]))
cube_shape304 = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength304))
spectra304 = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength304), dtype='float64', mode='r', shape=(cube_shape304[0], cube_shape304[1], cube_shape304[2]))
cube_shape171 = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength171))
spectra171 = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength171), dtype='float64', mode='r', shape=(cube_shape171[0], cube_shape171[1], cube_shape171[2]))
cube_shape193 = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength193))
spectra193 = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength193), dtype='float64', mode='r', shape=(cube_shape193[0], cube_shape193[1], cube_shape193[2]))
cube_shape211 = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength211))
spectra211 = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength211), dtype='float64', mode='r', shape=(cube_shape211[0], cube_shape211[1], cube_shape211[2]))

#global cpy_arr_spec
#cpy_arr_spec = np.copy(spectra)

global param1
param1 = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength1600))

"""
TIME = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))

if wavelength == 1600:
    time_step = 24
else:
    time_step = 12  # add as argument, or leave in as constant?
    
t_interp = np.linspace(0, TIME[len(TIME)-1], int(TIME[len(TIME)-1]/time_step))

if date == '20120923':   
    n_segments = 3  # break data into 12 segments of equal length
elif date == '20160426':
    n_segments = 12
elif date == '20140910':
    n_segments = 1
else:
    n_segments = 6
n = len(t_interp)
rem = n % n_segments  # n_segments is argument in module call (default to 1?)
freq_size = (n - rem) / n_segments
"""
### determine frequency values that FFT will evaluate

time_step = 24  
freq_size = (cube_shape1600[2]*2)+1
sample_freq24 = fftpack.fftfreq(freq_size, d=time_step)
pidxs24 = np.where(sample_freq24 > 0)    

time_step = 12
freq_size = (cube_shape304[2]*2)+1
sample_freq12 = fftpack.fftfreq(freq_size, d=time_step)
pidxs12 = np.where(sample_freq12 > 0)  

if 1:
    
    global f24
    global f12
    
    freqs24 = sample_freq24[pidxs24]
    freqs12 = sample_freq12[pidxs12]
    f24 = np.linspace(freqs24[0],freqs24[len(freqs24)-1],int(spectra1600.shape[2]))
    f12 = np.linspace(freqs12[0],freqs12[len(freqs12)-1],int(spectra304.shape[2]))
    print len(f24)
    print len(f12)
        
    global t1700, t1600, t304, t171, t193, t211
    t1700 = 0
    t1600 = 0
    t304 = 0
    t171 = 1
    t193 = 0
    t211 = 0
    
    
    h_map = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength1600))
 
    vis = np.load('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength1600))
    

    date_title = '%i/%02i/%02i' % (int(date[0:4]),int(date[4:6]),int(date[6:8]))

    

    # create list of titles and colorbar names for display on the figures
    titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Rolloever - [min]', 'Gaussian Amplitude', 'Gaussian Location -- [min]', 'Gaussian Width', 'F-Statistic', 'Visual Image - Averaged']
    
    # create figure with heatmap and spectra side-by-side subplots
    fig1 = plt.figure(figsize=(20,10))
    ax1 = plt.gca()
    ax1 = plt.subplot2grid((1,11),(0, 0), colspan=5, rowspan=1)
    plt.subplots_adjust(top=0.15)
    plt.subplots_adjust(left=0.25)
    ax1.set_xlim(0, h_map.shape[2]-1)
    ax1.set_ylim(0, h_map.shape[1]-1)  
    ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength1700, date_title, titles[1]), y = 1.01, fontsize=17)
    
    # was getting error "'AxesImage' object is not iterable"
    # - found: "Each element in img needs to be a sequence of artists, not a single artist."
    param = h_map[1]  # set initial heatmap to power law index     
    h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    im, = ([ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)])
    
    
    # design colorbar for heatmaps
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
    cbar.ax.tick_params(labelsize=13, pad=3)   
    
    
    # make toggle buttons to display each parameter's heatmap
    axcoeff = plt.axes([0.01, 0.9, 0.05, 0.063])
    axindex = plt.axes([0.07, 0.9, 0.05, 0.063])
    axroll = plt.axes([0.13, 0.9, 0.05, 0.063])
    axgauss_amp = plt.axes([0.19, 0.9, 0.05, 0.063])
    axgauss_loc = plt.axes([0.25, 0.9, 0.05, 0.063])
    axgauss_wid = plt.axes([0.31, 0.9, 0.05, 0.063])
    axfstat = plt.axes([0.37, 0.9, 0.05, 0.063])
    axvisual = plt.axes([0.43, 0.9, 0.05, 0.063])
    ax1700 = plt.axes([0.07, 0.05, 0.05, 0.063])
    ax1600 = plt.axes([0.13, 0.05, 0.05, 0.063])
    ax304 = plt.axes([0.19, 0.05, 0.05, 0.063])
    ax171 = plt.axes([0.25, 0.05, 0.05, 0.063])
    ax193 = plt.axes([0.31, 0.05, 0.05, 0.063])
    ax211 = plt.axes([0.37, 0.05, 0.05, 0.063])

    
 
    # set up spectra subplot
    ax2 = plt.subplot2grid((1,11),(0, 6), colspan=5, rowspan=1)
    ax2.loglog()
    ax2.set_xlim(10**-4.5, 10**-1.3)
    ax2.set_ylim(10**-5, 10**0)  
    
    fig1.canvas.mpl_connect('button_press_event', onclick)
    
    ax2.set_title('Spectra Fit', fontsize=15)
    plt.tight_layout()
    
    
    # add callbacks to each button - linking corresponding action
    callback = Index()
    
    bcoeff = Button(axcoeff, 'Coeff.')
    bcoeff.on_clicked(callback.coeff)
    bindex = Button(axindex, 'Index')
    bindex.on_clicked(callback.index)
    broll = Button(axroll, 'Rollover')
    broll.on_clicked(callback.roll)
    bgauss_amp = Button(axgauss_amp, 'Gauss Amp')
    bgauss_amp.on_clicked(callback.gauss_amp)
    bgauss_loc = Button(axgauss_loc, 'Gauss Loc')
    bgauss_loc.on_clicked(callback.gauss_loc)
    bgauss_wid = Button(axgauss_wid, 'Gauss Wid')
    bgauss_wid.on_clicked(callback.gauss_wid)
    bfstat = Button(axfstat, 'F-Stat')
    bfstat.on_clicked(callback.fstat)
    bvisual = Button(axvisual, 'Visual')
    bvisual.on_clicked(callback.visual)
    b1700 = Button(ax1700, '1700A')
    b1700.on_clicked(callback.f1700)
    b1600 = Button(ax1600, '1600A')
    b1600.on_clicked(callback.f1600)
    b304 = Button(ax304, '304A')
    b304.on_clicked(callback.f304)
    b171 = Button(ax171, '171A')
    b171.on_clicked(callback.f171)
    b193 = Button(ax193, '193A')
    b193.on_clicked(callback.f193)
    b211 = Button(ax211, '211A')
    b211.on_clicked(callback.f211)
    
    
plt.draw()