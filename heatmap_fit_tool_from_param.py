# -*- coding: utf-8 -*-
"""
Created on Mon Mar 06 07:34:15 2017

@author: Brendan
"""

#from __future__ import print_function  # not sure why this was included - keep in case
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
import sunpy
from sunpy.map import Map



def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)
    
    
class Index(object):
    ind = 0
    
       
    def pl_A(self, event):
        param = h_map[0]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max, picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[0]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[0], size=15, labelpad=10)
        plt.draw()

    def index(self, event):
        param = h_map[1]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
        plt.draw()
        
    def pl_C(self, event):
        param = h_map[2]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[2]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[2], size=15, labelpad=10)
        plt.draw()
        
    def gauss_amp(self, event):
        param = h_map[3]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[3]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[3], size=15, labelpad=10)
        plt.draw()
        
    def gauss_loc(self, event):
        param = 1./(np.exp(h_map[4]))
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[4]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[4], size=15, labelpad=10)
        plt.draw()
        
    def gauss_wid(self, event):
        param = h_map[5]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[5]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[5], size=15, labelpad=10)
        plt.draw()
        
    def f_stat(self, event):
        param = h_map[6]
        NaN_replace = np.nan_to_num(param)  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='none', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[6]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[6], size=15, labelpad=10)
        plt.draw()
        
    def visual(self, event):
        param = vis[0]
        im = ax1.imshow(param, cmap='sdoaia%i' % wavelength, interpolation='nearest', picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[7]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[7], size=15, labelpad=10)
        plt.draw()
        
    def roll(self, event):
        param = 1./((h_map[2]/h_map[0])**(-1./h_map[1]))
        im = ax1.imshow(param, cmap='jet_r', vmin=0, vmax=2000,interpolation='nearest', picker=True)
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[7]), y = 1.01, fontsize=17)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[7], size=15, labelpad=10)
        plt.draw()
        
    
    

# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy, c
    ix, iy = event.xdata, event.ydata
    ax2.clear()
    plt.draw()
    print ('x = %d, y = %d' % ( ix, iy))  # print location of pixel
    
    A2,n2,C2,P2,fp2,fw2 = h_map[0:6,iy,ix]

    m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)

    ax2.loglog(f,m2_fit)

    ax2.set_xlim(10**-4, 10**-1.3)
    ax2.set_ylim(10**-4.5, 10**0)  
    axvline(x=0.00333,color='k',ls='dashed', label='5 minutes')
    axvline(x=0.00555,color='k',ls='dotted', label='3 minutes')
    ax2.set_title('M2 Fit : Pixel (%ii , %ij)' % (iy, ix), fontsize=15)
    ax2.legend()
    plt.draw()
    
    return ix, iy
    


if 1:
    
    global f
    
    # determine frequency values that FFT will evaluate
    freq_size = ((299)*2) + 1  # determined from FFT-averaging script
    time_step = 12  # add as argument, or leave in as constant?
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    f = freqs  # frequencies
   
    
    global toggle
    toggle = 0
    
    
    vis = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_171_-500_500i_-500_600j_visual.npy')
    h_map = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_171_-500_500i_-500_600j_param_slope6_arthm.npy')
            
    """      
    # if I can make into a function -- use these, and take out corresponding lines below
    wavelength = wavelength
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    date_title = '%s-%s-%s' % (year,month,day)
    """
    
    #"""
    wavelength = 171
    year = 2013
    month = 6
    day = 26
    date_title = '%i-%02i-%02i' % (year,month,day)
    #"""   
    
    # arrays containing interesting points to be clicked for each dataset
    if wavelength == 211 and year == 2012 and month == 9 and day == 23:
        x = [250, 359, 567, 357, 322, 315, 97, 511, 316, 336]  
        y = [234, 308, 218, 197, 201, 199, 267, 5, 175, 181]
        
    if wavelength == 193 and year == 2013 and month == 5 and day == 30:
        x = [1]  
        y = [1]
    
    # create list of titles and colorbar names for display on the figures
    titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width', '$/chi^2$', 'Visual Image - Averaged']
    cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '$/chi^2$', 'Intensity']
    
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
    cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
    cbar.ax.tick_params(labelsize=13, pad=3)   
    
    
    # make toggle buttons to display each parameter's heatmap
    axpl_A = plt.axes([0.01, 0.9, 0.05, 0.075])
    axindex = plt.axes([0.07, 0.9, 0.05, 0.075])
    axpl_C = plt.axes([0.13, 0.9, 0.05, 0.075])
    axgauss_amp = plt.axes([0.19, 0.9, 0.05, 0.075])
    axgauss_loc = plt.axes([0.25, 0.9, 0.05, 0.075])
    axgauss_wid = plt.axes([0.31, 0.9, 0.05, 0.075])
    axf_stat = plt.axes([0.37, 0.9, 0.05, 0.075])
    axvisual = plt.axes([0.43, 0.9, 0.05, 0.075])
    axroll = plt.axes([0.49, 0.9, 0.05, 0.075])
 
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
    
    bpl_A = Button(axpl_A, 'pl_A')
    bpl_A.on_clicked(callback.pl_A)
    bindex = Button(axindex, 'Index')
    bindex.on_clicked(callback.index)
    bpl_C = Button(axpl_C, 'pl_C')
    bpl_C.on_clicked(callback.pl_C)
    bgauss_amp = Button(axgauss_amp, 'Gauss Amp')
    bgauss_amp.on_clicked(callback.gauss_amp)
    bgauss_loc = Button(axgauss_loc, 'Gauss Loc')
    bgauss_loc.on_clicked(callback.gauss_loc)
    bgauss_wid = Button(axgauss_wid, 'Gauss Wid')
    bgauss_wid.on_clicked(callback.gauss_wid)
    bf_stat = Button(axf_stat, 'F-Stat')
    bf_stat.on_clicked(callback.f_stat)
    bvisual = Button(axvisual, 'Visual')
    bvisual.on_clicked(callback.visual)
    broll = Button(axroll, 'Rollover')
    broll.on_clicked(callback.roll)
    
plt.draw()