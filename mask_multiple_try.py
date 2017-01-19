# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 09:36:34 2017

@author: Brendan
"""

# cant get param buttons and slider to work at same time

# you can pick param first, then move slider, but that is the only order where can use both

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from pylab import *
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import figure, show
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
with h5py.File("F:/Users/Brendan/Desktop/SolarProject/hdf5/20130530_193A_2300_2600i_2200_3000j_float_dogbox.hdf5",'r') as f:
   
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
            
        
    
    def update_min(val):
        vmin = s_vmin.val
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i][j] = h_map[1][i][j]
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                if R[i][j] < vmin:
                    R[i][j] = np.nan
     
        ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
        plt.subplots_adjust(top=0.85)
        plt.subplots_adjust(bottom=0.15)
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
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i][j] = h_map[1][i][j]
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                if R[i][j] > vmax:
                    R[i][j] = np.nan
    
        ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
        plt.subplots_adjust(top=0.85)
        plt.subplots_adjust(bottom=0.15)
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
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i][j] = h_map[1][i][j]
        ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)    

   
    if 1:
                    
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
        
            
        R = np.zeros((h_map.shape[1],h_map.shape[2]))



        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i][j] = h_map[1][i][j]
        
        # create list of titles and colorbar names for display on the figures
        titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width', '$/chi^2$', 'Visual Image - Averaged']
        cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '$/chi^2$', 'Intensity']
        
        # create figure with heatmap and spectra side-by-side subplots
        fig1 = plt.figure(figsize=(20,10))
        ax1 = plt.gca()
        ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
        plt.subplots_adjust(top=0.85)
        plt.subplots_adjust(bottom=0.15)
        ax1.set_xlim(0, h_map.shape[2]-1)
        ax1.set_ylim(0, h_map.shape[1]-1)  
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
        
        
        # was getting error "'AxesImage' object is not iterable"
        # - found: "Each element in img needs to be a sequence of artists, not a single artist."
        param = h_map[1]  # set initial heatmap to power law index     
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im, = ([ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)])
        
        vmin_initial = h_min
        vmax_initial = h_max
        
        
        # design colorbar for heatmaps
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
        cbar.ax.tick_params(labelsize=13, pad=3)   
        
        
        # make toggle buttons to display each parameter's heatmap
        axpl_A = plt.axes([0.27, 0.9, 0.05, 0.075])
        axindex = plt.axes([0.33, 0.9, 0.05, 0.075])
        axpl_C = plt.axes([0.39, 0.9, 0.05, 0.075])
        axgauss_amp = plt.axes([0.45, 0.9, 0.05, 0.075])
        axgauss_loc = plt.axes([0.51, 0.9, 0.05, 0.075])
        axgauss_wid = plt.axes([0.57, 0.9, 0.05, 0.075])
        axchi2 = plt.axes([0.63, 0.9, 0.05, 0.075])
        axvisual = plt.axes([0.69, 0.9, 0.05, 0.075])
        
        #axcolor = 'white'
        ax_vmin = plt.axes([0.25, 0.08, 0.45, 0.04])
        ax_vmax = plt.axes([0.25, 0.03, 0.45, 0.04])
        resetax = plt.axes([0.77, 0.055, 0.1, 0.04])
 
        
        #plt.tight_layout()
        
        
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
        s_vmin = Slider(ax_vmin, 'vmin', 0.1, 3., valinit=vmin_initial)
        s_vmax = Slider(ax_vmax, 'vmax', 0.1, 3., valinit=vmax_initial)
        s_vmin.on_changed(update_min)
        s_vmax.on_changed(update_max)
        button = Button(resetax, 'Reset')
        button.on_clicked(reset)
        
    plt.draw()

#plt.show()