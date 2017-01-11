# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 00:57:23 2017

@author: Brendan
"""

"""
############################
############################
# generate visual images
############################
############################
"""

## probably should include in heatmap generation function
## probably should have arrays included in param array - 
## created during 3x3 fitting?

## obviously need to make so that can load in different image arrays


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

def visual(dataset, date, wavelength, path_name):
    """
    Generates visual images for the region.
    
    dataset : 
        Path of file containing visual image data. (String)
        
    date : 
        Date of dataset in 'YYYYMMDD' format. (String)
    
    wavelength :
        Wavelength of dataset. (Integer)
                   
    path_name : 
        The directory to which the files should be saved. (String)
      
    Example:
    ::
        fm.visual(dataset='C:/Users/Brendan/Desktop/SDO/param_20120923_211_0_300_0_574_numpy.npy',
          wavelength=211, path_name='C:/Users/Brendan/Desktop/PHYS 326') 
    """

    titles = ['Average', 'Middle-File']
    names = ['average', 'mid']
    
    wavelength = wavelength
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    date_title = '%s-%s-%s' % (year,month,day)
    
    vis = np.load('%s' % dataset)
    
    for i in range(2):
        
        fig = plt.figure(figsize=(15,9))
        ax = plt.gca()
        plt.title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[i]), y = 1.01, fontsize=25)
        #im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
        im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength)
        plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
        plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('Intensity', size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=17, pad=5) 
        plt.tight_layout()
        plt.savefig('%s/%s_%i_visual_%s.jpeg' % (path_name, date, wavelength, names[i]))
        

def other(who, what, when):
    """
    Just a test example
    
    who : 
        Path of file containing visual image data. (String)
        
    what : 
        Date of dataset in 'YYYYMMDD' format. (String)
    
    when :
        Wavelength of dataset. (Integer)
                        
    Example:
    ::
        fm.visual(dataset='C:/Users/Brendan/Desktop/SDO/param_20120923_211_0_300_0_574_numpy.npy',
          wavelength=211, path_name='C:/Users/Brendan/Desktop/PHYS 326') 
    """

    titles = ['Average', 'Middle-File']
    names = ['average', 'mid']
    
    wavelength = wavelength
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    date_title = '%s-%s-%s' % (year,month,day)
