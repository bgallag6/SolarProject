# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 18:47:07 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import RectangleSelector
import matplotlib.patches as patches
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from matplotlib.widgets import Cursor
from pylab import axvline
import sunpy
from scipy import signal
from scipy import fftpack
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Button
from matplotlib.widgets import RadioButtons
from scipy.stats import f as ff
import os
from sunpy.map import Map

    
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def plotMap(p):
        param = h_map[p]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max, picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[p]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()

def plotMask(p):
        ax1.clear()
        ax1.set_xlim(0, h_map.shape[2]-1)
        ax1.set_ylim(0, h_map.shape[1]-1)  
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[p]), y = 1.01, fontsize=17)
        
        param = h_map[p]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        
        # generate p-value heatmap + masked Gaussian component heatmaps
        df1, df2 = 3, 6  # degrees of freedom for model M1, M2
        p_val = ff.sf(h_map[6], df1, df2)
        mask_thresh = 0.005  # significance threshold - masked above this value
        param_mask = np.copy(param) 
        param_mask[p_val > mask_thresh] = np.NaN  # mask the Gaussian component arrays with NaNs if above threshold
        im = ax1.imshow(param_mask, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max, picker=True)
        plt.colorbar(im,cax=cax)
        plt.draw()

    
class Index(object):
    ind = 0
         
    def coeff(self, event):
        global marker
        marker = 0
        plotMap(marker)  # could just do this
        return marker       

    def index(self, event):
        global marker
        marker = 1
        plotMap(marker)
        return marker     
        
    def roll(self, event):  # meh, should probably fix this
        global marker
        marker = 2
        paramA = h_map[0]
        paramn = h_map[1]
        paramC = h_map[2]
        param = (paramC/paramA)**(1./paramn)
        param = np.nan_to_num(param)/60.
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[2]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()      
        
    def gauss_amp(self, event):
        global marker
        marker = 3
        param = h_map[3]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[3]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
    
    def gauss_loc(self, event):
        global marker
        marker = 4
        param = (1./(np.exp(h_map[4]))/60.)
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet_r', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[4]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
        
    def gauss_wid(self, event):
        global marker
        marker = 5
        param = h_map[5]
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[5]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
        
    def fstat(self, event):
        global marker
        marker = 6
        param = h_map[6]
        NaN_replace = np.nan_to_num(param)  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        #h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        #h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        h_min = 0
        h_max = 30
        im = ax1.imshow(param, cmap='jet', interpolation='none', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[6]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
        
    def visual(self, event):
        param = vis[0]
        im = ax1.imshow(param, cmap='sdoaia%i' % wavelength, interpolation='nearest', picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[7]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
    
    def hist(self, event):
        global marker
        ax2.clear()
        ax2.set_title('Histogram: %s' % titles[marker], y = 1.01, fontsize=17)
        ax2.hist(np.reshape(h_map[marker], (h_map[marker].shape[0]*h_map[marker].shape[1])), bins=25)
        plt.draw()
    
    """    
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
    """
        
    def mask(self, event):
        global toggle2
        global marker
        if toggle2 == 0:
            toggle2 = 1
            plotMask(marker)
        elif toggle2 == 1:
            toggle2 = 0
            plotMap(marker)              
        return toggle2
        
    def saveFig(self, event):
        global count

        outdir = 'C:/Users/Brendan/Desktop/Tool_Figures/'

        if not os.path.exists(os.path.dirname(outdir)):
            try:
                print "Specified directory not found."
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST: raise
        else:
            plt.savefig('%s%s_%i_%i.pdf' % (outdir, date, wavelength, count), bbox_inches='tight')
        count += 1
        return count
  

# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy, c
    ixx, iyy = event.xdata, event.ydata
    if ixx > 1. and iyy > 1.:
        ax2.clear()
        del ax1.collections[:]
        plt.draw()
        print ('x = %d, y = %d' % ( ixx, iyy))  # print location of pixel
        ix = int(ixx)
        iy = int(iyy)
        
        s = np.zeros((spectra.shape[2]))
        m = np.zeros((spectra.shape[2]))
        g = np.zeros((spectra.shape[2]))
        pl = np.zeros((spectra.shape[2]))
    
        s[:] = spectra[iy][ix][:]
        
        m = GaussPowerBase(f_fit, param1[0][iy][ix], param1[1][iy][ix], param1[2][iy][ix], param1[3][iy][ix], param1[4][iy][ix], param1[5][iy][ix]) 
        g = Gauss(f_fit, param1[3][iy][ix], param1[4][iy][ix], param1[5][iy][ix])
        pl = PowerLaw(f_fit, param1[0][iy][ix], param1[1][iy][ix], param1[2][iy][ix])
        
        ax2.set_xlim(10**-4.5, 10**-1.3)
        ax2.set_ylim(10**-5, 10**0)
        
        ax2.loglog(f_fit, s, 'blue')
        ax2.loglog(f_fit, m, 'purple', label='M2 Combined')
        ax2.loglog(f_fit, g, 'g--', label='Gaussian')
        ax2.loglog(f_fit, pl, 'g', label='Power Law')
         
        axvline(x=0.00333,color='k',ls='dashed', label='5 minutes')
        axvline(x=0.00555,color='k',ls='dotted', label='3 minutes')
        ax2.set_title('Spectra Fit: Pixel (%ii , %ij)' % (iy, ix), y = 1.01, fontsize=17)
        plt.text(0.006, 10**-0.31, r'$A$ = {0:0.2e}'.format(h_map[0][iy][ix]), fontsize=23)
        plt.text(0.0061, 10**-0.51, r'$n$ = {0:0.2f}'.format(h_map[1][iy][ix]), fontsize=23)
        plt.text(0.006, 10**-0.73, r'$R$ = %0.1f [min]' % ((1./(h_map[2][iy][ix] / h_map[0][iy][ix])**(-1./ h_map[1][iy][ix]))/60.), fontsize=23)
        plt.text(0.0061, 10**-0.95, r'$\alpha$ = {0:0.2e}'.format(h_map[3][iy][ix]), fontsize=23)
        plt.text(0.0061, 10**-1.15, r'$\beta$ = {0:0.1f} [min]'.format((1./np.exp(h_map[4][iy][ix]))/60.), fontsize=23)
        plt.text(0.0061, 10**-1.35, r'$\sigma$ = {0:0.3f}'.format(h_map[5][iy][ix]), fontsize=23)
        #plt.legend(loc='lower left', prop={'size':20})
        legend = ax2.legend(loc='lower left', prop={'size':15}, labelspacing=0.35)
        ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
        for label in legend.get_lines():
                label.set_linewidth(2.0)  # the legend line width
        plt.draw()

    return ix, iy
    
# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    #return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)    
    return A2*f2**-n2 + C2 + P2*(1./ ((np.pi*fw2)*(1.+((np.log(f2)-fp2)/fw2)**2)))

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
        
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    #return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2) 
    return P*(1./ ((np.pi*fw)*(1.+((np.log(f)-fp)/fw)**2)))
    

"""
##############################################################################
##############################################################################
"""

directory = 'F:'
#date = '20160327'
date = '20130626'
wavelength = 171

global spectra
global param1

cube_shape = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength))
spectra = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))

param1 = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))


global marker
global toggle
global toggle2
global count
marker = 1
toggle = 0
toggle2 = 0
count = 0


### determine frequency values that FFT will evaluate
if wavelength == 1600 or wavelength == 1700:
    time_step = 24
else:
    time_step = 12
freq_size = (cube_shape[2]*2)+1
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)    


if 1:
    
    global f_fit
    
    freqs = sample_freq[pidxs]
    print len(freqs)
    f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],int(spectra.shape[2]))   
    
    
    h_map = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
 
    vis = np.load('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength))

    
    # arrays containing interesting points to be clicked for each dataset
    if date == '20120923' and wavelength == 211:
        x = [250, 359, 567, 357, 322, 315, 97, 511, 316, 336]  
        y = [234, 308, 218, 197, 201, 199, 267, 5, 175, 181]
    
    # create list of titles and colorbar names for display on the figures
    titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Rollover [min]', 'Gaussian Amplitude', 'Gaussian Location [min]', 'Gaussian Width', 'F-Statistic', 'Averaged Visual Image']
    date_title = '%i/%02i/%02i' % (int(date[0:4]),int(date[4:6]),int(date[6:8]))
    
    # create figure with heatmap and spectra side-by-side subplots
    fig1 = plt.figure(figsize=(20,10))
    ax1 = plt.gca()
    ax1 = plt.subplot2grid((30,31),(4, 1), colspan=14, rowspan=25)

    ax1.set_xlim(0, h_map.shape[2]-1)
    ax1.set_ylim(0, h_map.shape[1]-1)  
    ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[1]), y = 1.01, fontsize=17)
    
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
    #axscatter = plt.axes([0.49, 0.9, 0.05, 0.063])
    axhist = plt.axes([0.49, 0.9, 0.05, 0.063])
    axmask = plt.axes([0.55, 0.9, 0.05, 0.063])
    axsaveFig = plt.axes([0.61, 0.9, 0.05, 0.063])
 
    # set up spectra subplot
    ax2 = plt.subplot2grid((30,31),(4, 16), colspan=14, rowspan=25)
    ax2.loglog()
    ax2.set_xlim(10**-4.5, 10**-1.3)
    ax2.set_ylim(10**-5, 10**0)  
    
    fig1.canvas.mpl_connect('button_press_event', onclick)
    
    ax2.set_title('Spectra Fit', y = 1.01, fontsize=17)
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
    #bscatter = Button(axscatter, 'Scatter')
    #bscatter.on_clicked(callback.scatter)
    bhist = Button(axhist, 'Hist.')
    bhist.on_clicked(callback.hist)
    bmask = Button(axmask, 'Mask')
    bmask.on_clicked(callback.mask)
    bsaveFig = Button(axsaveFig, 'Save')
    bsaveFig.on_clicked(callback.saveFig)
    
plt.draw()

