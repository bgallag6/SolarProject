# -*- coding: utf-8 -*-
"""
Created on Thu Feb 08 15:48:53 2018

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from sunpy.map import Map
import matplotlib


directory = 'D:'
#directory = 'C:/Users/Brendan/Desktop'
date = '20130626'
wavelengths = [171,193,304,1700]
letters = ['d','c','e','f']

tickmin = [20,20,10,660]
tickincr = [45,70,5,100]

matplotlib.rc('text', usetex = True)  # use with latex commands
#plt.rc('font', family='serif')
plt.rc('font',**{'family':'serif','serif':['Times']})

# in case font weight appears bold
#del matplotlib.font_manager.weight_dict['roman']
#matplotlib.font_manager._rebuild()

for q in range(4):
    wavelength = wavelengths[q]
    
    # load parameter array and visual images from file tree structure 
    heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
    visual = np.load('%s/DATA/Output/%s/%i/visual.npy'% (directory, date, wavelength))  
    
    visual = visual[1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
    h_map = heatmaps    
    
    #plt.rcParams["font.family"] = "Times New Roman"
    font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels
    
    wavelength = wavelength    
    
    # trim x/y dimensions equally so that resulting region is 1600x1600    
    trim_y = int((h_map.shape[1]-1600)/2)
    trim_x = int((h_map.shape[2]-1600)/2)
    h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)
    #h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x*2:]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    
    
    x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
    y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
    x_ind = [-800,-600,-400,-200,0,200,400,600,800]
    y_ind = [800,600,400,200,0,-200,-400,-600,-800]    
    
    """
    # for creating sunspot umbra + PPV contour overlays from 1600
    v1600 = np.load('%s/DATA/Output/%s/1600/visual.npy'% (directory, date))  
    v1600 = v1600[:,1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
    p1600 = np.load('%s/DATA/Output/%s/1600/param.npy'% (directory, date)) 
    
    h_map = h_map[:,:p1600.shape[1],:p1600.shape[2]]
    visual = visual[:v1600.shape[1],:v1600.shape[2]]
    """    
    
    fig_width = 10+2  # works better for 20130626 (with no x/y labels)
    fig_height = 10  # works better for 20130626
    
      
    # generate visual images
    titles_vis = ['Average', 'Middle-File']
    names_vis = ['average', 'mid']
    
    vis = visual
    
    trim_yv = int((vis.shape[0]-1600)/2)
    trim_xv = int((vis.shape[1]-1600)/2)
    vis = vis[trim_yv:vis.shape[0]-trim_yv, trim_xv:vis.shape[1]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides) 
    #vis = vis[trim_yv:vis.shape[0]-trim_yv, trim_xv*2:]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides) 
           
        
    v_min = np.percentile(vis,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    v_max = np.percentile(vis,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)  

    
    h_range = np.abs(v_max-v_min)
    h_step = h_range / 10.
    c_ticks = np.zeros((11))
    for h in range(11):
        c_ticks[h] = v_min + h_step*h 
    

    ## Panel A for each wavelength
    fig = plt.figure(figsize=(fig_width,fig_height))
    
    ax = plt.gca()
    plt.title('(a) Visual %s' % (titles_vis[0]), y = 1.02, fontsize=font_size)  # no date / wavelength   
    im = ax.imshow(np.flipud(vis), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
    plt.xticks(x_ticks,x_ind,fontsize=font_size)
    plt.yticks(y_ticks,y_ind,fontsize=font_size)
    ax.tick_params(axis='both', which='major', pad=10)
    
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax,format='%d')
    cbar.set_label('DN/s', size=font_size, labelpad=10)

    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)

    c_ind = np.array([tickmin[q]+(tickincr[q]*h) for h in range(11)])
    cbar.set_ticklabels(c_ind)
    
    #plt.savefig('C:/Users/Brendan/Desktop/%s_%i_visual_normF.pdf' % (date, wavelength), format='pdf', bbox_inches='tight')
    
      
    
    ## 6-panel visual average summary
    fig = plt.figure(figsize=(fig_width,fig_height))
    
    ax = plt.gca()
    plt.title('(%s) %i \AA' % (letters[q], wavelength), y = 1.02, fontsize=font_size)  # no date / wavelength
    im = ax.imshow(np.flipud(vis), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
    plt.xticks(x_ticks,x_ind,fontsize=font_size)
    plt.yticks(y_ticks,y_ind,fontsize=font_size)
    ax.tick_params(axis='both', which='major', pad=10)
    
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax,format='%d')
    cbar.set_label('DN/s', size=font_size, labelpad=10)

    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)
    
    c_ind = np.array([tickmin[q]+(tickincr[q]*h) for h in range(11)])
    cbar.set_ticklabels(c_ind)
    
    plt.savefig('C:/Users/Brendan/Desktop/%s_%i_visual_normF2.pdf' % (date, wavelength), format='pdf', bbox_inches='tight')