# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 15:02:12 2017

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



#def heatmap_spectra_tool(dataset, date, wavelength, num_segments):
#with h5py.File("%s" % dataset,'r') as f:    
with h5py.File("F:/Users/Brendan/Desktop/SolarProject/20130815_193_1000_1600i_1950_2950j_rebin2_FULL.hdf5",'r') as f:
   
    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]
        
        
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
            param = h_map[4]
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
            
        def chi2(self, event):
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
        m = np.zeros((m2_fit.shape[2]))
        p = np.zeros((powerlaw.shape[2]))
        g = np.zeros((gaussian.shape[2]))
        for i in range(0,cpy_arr_spec.shape[2]):
            s[i] = cpy_arr_spec[iy][ix][i]
        for i in range(0,cpy_arr_m2.shape[2]):
            m[i] = cpy_arr_m2[iy][ix][i]
        for i in range(0,cpy_arr_powerlaw.shape[2]):
            p[i] = cpy_arr_powerlaw[iy][ix][i]
        for i in range(0,cpy_arr_gaussian.shape[2]):
            g[i] = cpy_arr_gaussian[iy][ix][i]
        ax2.loglog(f_fit_full,s)
        ax2.loglog(f_fit,m)
        ax2.loglog(f_fit,p)
        ax2.loglog(f_fit,g)

        ax2.set_xlim(10**-4.5, 10**-1.3)
        ax2.set_ylim(10**-5, 10**0)  
        axvline(x=0.00333,color='k',ls='dashed', label='5 minutes')
        axvline(x=0.00555,color='k',ls='dotted', label='3 minutes')
        ax2.set_title('Spectra Fit : Pixel (%ii , %ij)' % (iy, ix), fontsize=15)
        plt.legend(loc='upper right', prop={'size':15})
        plt.draw()
        
        return ix, iy
        



    #m2 = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/M2_20141025_304_2300_2600_1700_2400_numpy.npy')
    #m2_fit = f['m2_fit']
    m2_fit = np.array(f['m2_fit'])
    powerlaw = np.array(f['powerlaw'])
    gaussian = np.array(f['gaussian'])

    
    #spectra = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/spectra_20141025_304_2300_2600_1700_2400_numpy.npy')
    #spectra = f['spectra']
    spectra = np.array(f['spectra'])

    global cpy_arr_spec
    global cpy_arr_m2
    global cpy_arr_powerlaw
    global cpy_arr_gaussian
    cpy_arr_spec = np.copy(spectra)
    cpy_arr_m2 = np.copy(m2_fit)
    cpy_arr_powerlaw = np.copy(powerlaw)
    cpy_arr_gaussian = np.copy(gaussian)
    
    TIME = f['time']  # need for now, when have m2 be 1/2 of num freqs, could just do m2.shape x 2 then
    # or just have argument of 3,6,12 hours, and either number of segments or hours per segment
    #TIME = np.load('F:/SDO/datacubes_numpy/SDO_20120923_211A_(528)_(132)x_(100)_100y_time.npy')
    
    t_interp = np.linspace(0, TIME[len(TIME)-1], TIME[len(TIME)-1]/12)

       
    n_segments = 3  # break data into 12 segments of equal length
    n = len(t_interp)
    rem = n % n_segments  # n_segments is argument in module call (default to 1?)
    freq_size = (n - rem) / n_segments
    
    ### determine frequency values that FFT will evaluate
    time_step = 12  # add as argument, or leave in as constant?
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)    
    
    
    if 1:
        
        global f_fit
        global f_fit_full
        
        freqs = sample_freq[pidxs]
        
        f_fit_full = np.linspace(freqs[0],freqs[len(freqs)-1],spectra.shape[2])
        f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],m2_fit.shape[2])
        
        
        global toggle
        toggle = 0
        

        #Y = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20141025_304_2300_2600_1700_2400_numpy.npy')
        #Y = f['params']  # reason why called within loop?
        h_map = np.array(f['params'])
        print(h_map.shape)
        
        vis = np.array(f['visual'])
        
                
        """      
        # if I can make into a function -- use these, and take out corresponding lines below
        wavelength = wavelength
        year = date[0:4]
        month = date[4:6]
        day = date[6:8]
        date_title = '%s-%s-%s' % (year,month,day)
        """
        
        #"""
        wavelength = 211
        year = 2012
        month = 9
        day = 23
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
        axchi2 = plt.axes([0.37, 0.9, 0.05, 0.075])
        axvisual = plt.axes([0.43, 0.9, 0.05, 0.075])
        axscatter = plt.axes([0.49, 0.9, 0.05, 0.075])
 
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
        bchi2 = Button(axchi2, 'Chi2')
        bchi2.on_clicked(callback.chi2)
        bvisual = Button(axvisual, 'Visual')
        bvisual.on_clicked(callback.visual)
        bscatter = Button(axscatter, 'Scatter')
        bscatter.on_clicked(callback.scatter)
        
    plt.draw()