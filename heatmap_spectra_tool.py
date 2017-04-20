# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 00:21:29 2017

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
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max, picker=True)
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
        param = h_map[6]
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
        
    def scatter(self, event):
        global toggle
        global col
        if toggle == 0:
            toggle = 1
            col = ax1.scatter(x, y, s=50, c='white', picker=True)
        elif toggle == 1:
            toggle = 0
            col.remove()
            plt.draw()                
        return toggle
        
  

# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy, c
    ix, iy = event.xdata, event.ydata
    ax2.clear()
    plt.draw()
    print ('x = %d, y = %d' % ( ix, iy))  # print location of pixel
    
    s = np.zeros((spectra.shape[2]))
    m = np.zeros((spectra.shape[2]))
    g = np.zeros((spectra.shape[2]))
    pl = np.zeros((spectra.shape[2]))
    
    for i in range(0,cpy_arr_spec.shape[2]):
        s[i] = cpy_arr_spec[iy][ix][i]

    m = GaussPowerBase(f_fit, param1[0][iy][ix], param1[1][iy][ix], param1[2][iy][ix], param1[3][iy][ix], param1[4][iy][ix], param1[5][iy][ix]) 
    g = Gauss(f_fit, param1[3][iy][ix], param1[4][iy][ix], param1[5][iy][ix])
    pl = PowerLaw(f_fit, param1[0][iy][ix], param1[1][iy][ix], param1[2][iy][ix])
    
    ax2.loglog(f_fit, s, 'blue')
    ax2.loglog(f_fit, m, 'purple', label='M2 Combined')
    ax2.loglog(f_fit, g, 'g--', label='Gaussian')
    ax2.loglog(f_fit, pl, 'g', label='Power Law')
     
    ax2.set_xlim(10**-4.5, 10**-1.3)
    ax2.set_ylim(10**-5, 10**0)  
    axvline(x=0.00333,color='k',ls='dashed', label='5 minutes')
    axvline(x=0.00555,color='k',ls='dotted', label='3 minutes')
    ax2.set_title('Spectra Fit : Pixel (%ii , %ij)' % (iy, ix), fontsize=15)
    plt.text(0.006, 10**-0.31, r'$A$ = {0:0.2e}'.format(h_map[0][iy][ix]), fontsize=23)
    plt.text(0.0061, 10**-0.51, r'$n$ = {0:0.2f}'.format(h_map[1][iy][ix]), fontsize=23)
    plt.text(0.006, 10**-0.73, r'$R$ = %0.1f [min]' % ((1./(h_map[2][iy][ix] / h_map[0][iy][ix])**(-1./ h_map[1][iy][ix]))/60.), fontsize=23)
    plt.text(0.0061, 10**-0.95, r'$\alpha$ = {0:0.2e}'.format(h_map[3][iy][ix]), fontsize=23)
    plt.text(0.0061, 10**-1.15, r'$\beta$ = {0:0.1f} [min]'.format((1./np.exp(h_map[4][iy][ix]))/60.), fontsize=23)
    plt.text(0.0061, 10**-1.35, r'$\sigma$ = {0:0.3f}'.format(h_map[5][iy][ix]), fontsize=23)
    #plt.legend(loc='lower left', prop={'size':20})
    legend = ax2.legend(loc='lower left', prop={'size':15}, labelspacing=0.35)
    for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
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

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20130626'
wavelength = 1600


cube_shape = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength))
spectra = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))

global cpy_arr_spec
cpy_arr_spec = np.copy(spectra)

global param1
param1 = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))

TIME = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))

t_interp = np.linspace(0, TIME[len(TIME)-1], TIME[len(TIME)-1]/12)

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

### determine frequency values that FFT will evaluate
if wavelength == 1600:
    time_step = 24
else:
    time_step = 12  # add as argument, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)    


if 1:
    
    global f_fit
    
    freqs = sample_freq[pidxs]
    
    f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],spectra.shape[2])
        
    global toggle
    toggle = 0
    
    
    h_map = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
 
    vis = np.load('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength))
    

    date_title = '%i/%02i/%02i' % (int(date[0:4]),int(date[4:6]),int(date[6:8]))

    
    # arrays containing interesting points to be clicked for each dataset
    if date == '20120923' and wavelength == 211:
        x = [250, 359, 567, 357, 322, 315, 97, 511, 316, 336]  
        y = [234, 308, 218, 197, 201, 199, 267, 5, 175, 181]
        
    if date == '20130530' and wavelength == 193:
        x = [1]  
        y = [1]
    
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
    ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
    
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
    axscatter = plt.axes([0.49, 0.9, 0.05, 0.063])
 
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
    bscatter = Button(axscatter, 'Scatter')
    bscatter.on_clicked(callback.scatter)
    
plt.draw()