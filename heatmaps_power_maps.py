# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:43:04 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm

def heatmap(directory, date, wavelength):
    
    # load parameter array and visual images from file tree structure 
    heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
    h_map = heatmaps
    
    plt.rcParams["font.family"] = "Times New Roman"
    font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels
    
    wavelength = wavelength    
    
    #"""
    # *for 2013/06/26 region* -- trim x/y dimensions equally so that resulting region is 1600x1600
    trim_y = (h_map.shape[1]-1600)/2
    trim_x = (h_map.shape[2]-1600)/2
    h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    
    
    x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
    y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
    x_ind = [-800,-600,-400,-200,0,200,400,600,800]
    y_ind = [800,600,400,200,0,-200,-400,-600,-800]    
    #"""
    
    """
    xdim = int(np.floor(h_map.shape[2]/100))
    ydim = int(np.floor(h_map.shape[1]/100))
    
    x_ticks = [100*i for i in range(xdim+1)]
    y_ticks = [100*i for i in range(ydim+1)]
    
    x_ind = x_ticks
    y_ind = y_ticks
    """    
    
    # generate p-value heatmap + masked Gaussian component heatmaps
    df1, df2 = 3, 6  # degrees of freedom for model M1, M2
    p_val = ff.sf(h_map[6], df1, df2)
    
    mask_thresh = 0.005  # significance threshold - masked above this value
      
    loc_mask = np.copy(h_map[4])

    # mask the Gaussian component arrays with NaNs if above threshold 
    loc_mask[p_val > mask_thresh] = np.NaN # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
    
    # determine percentage of region masked                     
    loc_mask = (1./np.exp(loc_mask))/60.  # convert Gaussian location to minutes
    count = np.count_nonzero(np.isnan(loc_mask))   
    total_pix = loc_mask.shape[0]*loc_mask.shape[1]
    mask_percent0 = ((np.float(count))/total_pix)*100
    
    for r in range(6):
        gmin = 3.6+(0.2*r)
        gmax = 3.8+(0.2*r)
        loc_41 = np.copy(loc_mask)
        loc_41 = np.reshape(np.copy(loc_mask),(loc_41.shape[0]*loc_41.shape[1]))
        loc_41 = [np.NaN if x < gmin else x for x in loc_41]
        loc_41 = [np.NaN if x > gmax else x for x in loc_41]
        loc_41 = np.reshape(loc_41, (loc_mask.shape[0], loc_mask.shape[1]))
        count = np.count_nonzero(np.isnan(loc_41))   
        total_pix = loc_41.shape[0]*loc_41.shape[1]
        mask_percent = (1-((np.float(count))/total_pix))*100
        plt.figure(figsize=(12,10))
        ax = plt.gca()
        plt.title(r'Gauss. Loc. | %0.1f-%0.1f min | $f_{coverage}$ = %0.1f%s' % (gmin,gmax,mask_percent,'%'), y=1.02, fontsize=font_size, fontname="Times New Roman")
        cmap = cm.get_cmap('Greys', 1)
        plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
        plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
        ax.tick_params(axis='both', which='major', pad=10)
        plt.imshow(np.flipud(loc_41), cmap=cmap)
        #plt.savefig('C:/Users/Brendan/Desktop/powermap_%imin.pdf' % (gmin*10))
    
    """
    ## works for other, non-square regions
    if h_map.shape[2] > h_map.shape[1]:
        aspect_ratio = float(h_map.shape[2]) / float(h_map.shape[1])
        fig_height = 10
        fig_width = 10*aspect_ratio
        
    else:
        aspect_ratio = float(h_map.shape[1]) / float(h_map.shape[2])
        #print aspect_ratio
        #fig_width = 10
        fig_width = 10+2  # works better for 20130626 (with no x/y labels)
        fig_height = 10*aspect_ratio  # works better for 20130626
    """
    
    #"""
    ## works for 2013/06/26 -- or square regions
    fig_width = 10+2  # works better for 20130626 (with no x/y labels)
    #fig_width = 10+3  # works better for 20130626 (with x/y labels)
    fig_height = 10  # works better for 20130626
    #"""
    
    i = 4  # index of gaussian location in parameter array
    
    ## Gaussian Location full heatmap
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    plt.title(r'Gauss. Loc. $\beta$ [min]', y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength

    h_map[i] = (1./(np.exp(h_map[i])))/60.
    h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    cmap = cm.get_cmap('jet_r', 10)

    # specify colorbar ticks to be at boundaries of segments
    h_range = np.abs(h_max-h_min)
    h_step = h_range / 10.
    c_ticks = np.zeros((11))
    for h in range(11):
        c_ticks[h] = h_min + h_step*h 
        
    im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=h_min, vmax=h_max)
    plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
    plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax, format='%0.1f')
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.jpeg' % (directory, date, wavelength, date, wavelength, names[i]))
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf')
    
    
    ## Gaussian Location masked heatmap
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar
    plt.title(r'Gauss. Loc. [min] $\beta$; $p$ < %0.3f | $f_{masked}$ = %0.1f%s' % (mask_thresh, mask_percent0,'%'), y = 1.02, fontsize=font_size, fontname="Times New Roman")
    cmap = cm.get_cmap('jet_r', 10)
                  
    im = ax.imshow(np.flipud(loc_mask), cmap = cmap, vmin=h_min, vmax=h_max)
    plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
    plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax, format='%0.1f')
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s_mask_%i.jpeg' % (directory, date, wavelength, date, wavelength, names[i], (1./mask_thresh)))
    #plt.savefig('C:/Users/Brendan/Desktop/%s_%i_%s_mask_%i.pdf' % (date, wavelength, names[i], (1./mask_thresh)), format='pdf')
    
    
    ## Gaussian Location histogram
    flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
    
    # calculate some statistics
    mean = np.mean(flat_param)
    sigma = np.std(flat_param)   
    
    fig = plt.figure(figsize=(fig_width+1,fig_height))
    plt.title(r'Gauss. Loc. $\beta$', y = 1.02, fontsize=font_size)  # no date / wavelength
    plt.xlabel('Period [min]', fontsize=font_size, labelpad=10)
    plt.ylabel('Bin Count', fontsize=font_size, labelpad=10)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim(h_min, h_max)
    y, x, _ = plt.hist(flat_param, bins=200, range=(h_min, h_max))
    #n, bins, patches = plt.hist(flat_param, bins=200, range=(h_min, h_max))
    n=y[1:-2]
    bins=x[1:-2]
    elem = np.argmax(n)
    bin_max = bins[elem]
    plt.vlines(bin_max, 0, y.max()*1.1, color='blue', linestyle='dotted', linewidth=1.5, label='mode=%0.6f' % bin_max)      
    #y, x, _ = plt.hist(flat_param, bins=100)
    plt.ylim(0, y.max()*1.1)
    plt.vlines(mean, 0, y.max()*1.1, color='red', linestyle='solid', linewidth=1.5, label='mean=%0.6f' % mean)
    plt.vlines(0, 0, y.max()*1.1, color='white', linestyle='dashed', linewidth=1.5, label='sigma=%0.6f' % sigma) 
    #plt.vlines(mode[0], 0, y.max()*1.1, color='red', linestyle='dashed', linewidth=1.5, label='mode=%0.6f' % mode[0])
    #plt.vlines(lower_sigma, 0, y.max()*1.1, color='blue', linestyle='dotted', linewidth=1.5, label='low 68=%0.6f' % lower_sigma)
    #plt.vlines(upper_sigma, 0, y.max()*1.1, color='blue', linestyle='dashed', linewidth=1.5, label='high 68=%0.6f' % upper_sigma)
    legend = plt.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
    for label in legend.get_lines():
        label.set_linewidth(2.0)  # the legend line width
    
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_Histogram_%s.jpeg' % (directory, date, wavelength, date, wavelength, names[i]))
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_Histogram_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf')
    #"""



directory = 'D:/Users/Brendan/Desktop/SolarProject'
date = '20130626'
wavelength = 1600

heatmap(directory= '%s' % (directory), date='%s' % (date), wavelength= wavelength)