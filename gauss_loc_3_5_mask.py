# -*- coding: utf-8 -*-
"""
Created on Sat Apr 01 21:35:14 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = 20140818
wavelength = 1600
    
# create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'Power Law Index - $n$', r'Power Law Tail - $C$', r'Gaussian Amplitude [flux] - $\alpha$', r'Gaussian Location [min] - $\beta$', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', 'p-Value']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [min]', 'Width', 'F-Statistic', 'Amplitude Scaled', 'p-Value']

#vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore  (or option to set ranges for specific wavelengths?)
#vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore

# load parameter array and visual images from file tree structure 
heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
visual = np.load('%s/DATA/Output/%s/%i/visual.npy'% (directory, date, wavelength))  
h_map = heatmaps

#heatmaps = np.load('%s/DATA/Output/%s/%i/*param.npy' % (directory, date, wavelength))
#visual = np.load('%s/DATA/Output/%s/%i/*visual.npy'% (directory, date, wavelength))    

font_size = 23  # set the font size to be used for all text - titles, tick marks, text, labels

wavelength = wavelength    

#h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]  # trim last row and column from array (originally needed since went one past)

# trim x/y dimensions equally so that resulting region is 1600x1600    
trim_y = (h_map.shape[1]-1600)/2
trim_x = (h_map.shape[2]-1600)/2
#h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

#h_map = h_map[:, 0:h_map.shape[1]-50, 0:500]  # for 20130626 blobs      
#x_ticks = [0,100,200,300,400,500]
#y_ticks = [0,100,200,300,400]   
#y_ticks = [0,100,200,300] # 20120923

#x_ticks = [0,100,200,300,400,500]
#y_ticks = [0,100,200,300]  

x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
x_ind = [-800,-600,-400,-200,0,200,400,600,800]
y_ind = [800,600,400,200,0,-200,-400,-600,-800]


# generate p-value heatmap + masked Gaussian component heatmaps
df1, df2 = 3, 6  # degrees of freedom for model M1, M2
p_val = ff.sf(h_map[6], df1, df2)

p_mask = np.copy(p_val)

mask_thresh = 0.005  # significance threshold - masked above this value
   
p_mask = np.copy(p_val)
amp_mask = np.copy(h_map[3])
loc_mask = np.copy(h_map[4])
wid_mask = np.copy(h_map[5])

# mask the Gaussian component arrays with NaNs if above threshold 
p_mask[p_val > mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
amp_mask[p_val > mask_thresh] = np.NaN
loc_mask[p_val > mask_thresh] = np.NaN
wid_mask[p_val > mask_thresh] = np.NaN

# determine percentage of region masked 
count = np.count_nonzero(np.isnan(p_mask))   
total_pix = p_val.shape[0]*p_val.shape[1]
mask_percent = ((np.float(count))/total_pix)*100
                
loc_mask = (1./np.exp(loc_mask))/60.  # convert Gaussian location to minutes

loc_mask[loc_mask > 6.2] = np.NaN
loc_mask[loc_mask < 2.0] = np.NaN

plots = [p_mask, amp_mask, loc_mask, wid_mask]  # make array of masked plots to iterate over


#"""
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
#"""

#fig_width = 10+2  # works better for 20130626 (with no x/y labels)
#fig_width = 10+3  # works better for 20130626 (with x/y labels)
#fig_height = 10  # works better for 20130626

#for i in range(0,len(titles)-1):
for i in range(4,5):
    
    #fig = plt.figure(figsize=(13,9))
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
    plt.title('%s' % (titles[i]), y = 1.02, fontsize=font_size)  # no date / wavelength
    
    if i == 6:
        NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        #cmap = 'jet'
        cmap = cm.get_cmap('jet', 10)  # specify discrete colorscale with 10 intervals 
    elif i == 4:
        h_map[i] = (1./(np.exp(h_map[i])))/60.
        h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        #cmap = 'jet_r'  # reverse color-scale for Gaussian Location, because of flipped frequencies to seconds
        cmap = cm.get_cmap('jet_r', 10)
    elif i == 8:
        df1, df2 = 3, 6
        h_map[6] = ff.sf(h_map[6], df1, df2)
        h_min = np.percentile(h_map[6],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(h_map[6],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        #cmap = 'jet'      
        cmap = cm.get_cmap('jet', 10)               
    else:
        h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        #cmap = 'jet'
        cmap = cm.get_cmap('jet', 10)
    
    # specify colorbar ticks to be at boundaries of segments
    h_range = np.abs(h_max-h_min)
    h_step = h_range / 10.
    #h1 = h_min + (h_step/2.)
    c_ticks = np.zeros((11))  # use 10 if centers of segments
    for h in range(11):
        #c_ticks[h] = h1 + h_step*h  # use for centers of segments
        c_ticks[h] = h_min + h_step*h 
        
    im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=h_min, vmax=h_max)
    #im = ax.imshow(h_map[i], cmap = cmap, vmin=h_min, vmax=h_max)
    #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
    #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
    #plt.xticks(x_ticks,fontsize=font_size)
    #plt.yticks(y_ticks,fontsize=font_size)
    #plt.xticks(x_ticks,x_ind,fontsize=font_size)
    #plt.yticks(y_ticks,y_ind,fontsize=font_size)
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    if i == 4:
        cbar = plt.colorbar(im,cax=cax, format='%0.1f')
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    #cbar.set_ticks(np.round(c_ticks,8))  # 8 for slope (or might as well be for all, format separately)
    cbar.set_ticks(c_ticks)  # 8 for slope (or might as well be for all, format separately)
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.jpeg' % (directory, date, wavelength, date, wavelength, names[i]))
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf')
    
    num_grad = 15    
    
    if i == 2 or i == 3 or i == 4 or i == 5:   
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        if i == 2:
            plt.title(r'$p-value$ < %0.3f | f_{masked} = %0.1f' % (mask_thresh, mask_percent), y = 1.02, fontsize=font_size)
        else:
            plt.title(r'%s: $p$ < %0.3f | f_{masked} = %0.1f' % (titles[i], mask_thresh, mask_percent), y = 1.02, fontsize=font_size)
        if i == 4:
            cmap = cm.get_cmap('jet_r',num_grad)
            par = np.nan_to_num(plots[i-2])
            h_min = 2.  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = 6.2  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            h_range = np.abs(h_max-h_min)
            h_step = h_range / num_grad
            #h1 = h_min + (h_step/2.)
            c_ticks = np.zeros((num_grad+1))
            for h in range(num_grad+1):
                #c_ticks[h] = h1 + h_step*h
                c_ticks[h] = h_min + h_step*h
                
            
        im = ax.imshow(np.flipud(plots[i-2]), cmap = cmap, vmin=h_min, vmax=h_max)
    
        #h_step = h_range / 10.
        #h1 = h_min + (h_step/2.)
        #c_ticks = np.zeros((10))
        #for h in range(10):
        #    c_ticks[h] = h1 + h_step*h
        #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
        #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
        #plt.xticks(x_ticks,fontsize=font_size)
        #plt.yticks(y_ticks,fontsize=font_size)
        #plt.xticks(x_ticks,x_ind,fontsize=font_size)
        #plt.yticks(y_ticks,y_ind,fontsize=font_size)
        ax.tick_params(axis='both', which='major', pad=10)
        divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
        cax = divider.append_axes("right", size="3%", pad=0.07)
        if i == 2:
            cbar = plt.colorbar(im,cax=cax, format='%0.4f')
        elif i == 3:
            cbar = plt.colorbar(im,cax=cax, format='%0.4f')
        elif i == 4:
            cbar = plt.colorbar(im,cax=cax, format='%0.1f')
        elif i == 5:
            cbar = plt.colorbar(im,cax=cax, format='%0.2f')
        #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=font_size, pad=5) 
        #cbar.set_ticks(np.round(c_ticks,8))  # 8 for slope
        cbar.set_ticks(c_ticks)  # 8 for slope
        #plt.tight_layout()
        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_m[k]))
        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s_mask_%i.pdf' % (directory, date, wavelength, date, wavelength, names[i], (1./mask_thresh)), format='pdf')
    
    #"""
    flat_param = np.reshape(plots[i-2], (plots[i-2].shape[0]*plots[i-2].shape[1]))

    fig = plt.figure(figsize=(fig_width+1,fig_height))
    plt.title('%s' % (titles[i]), y = 1.02, fontsize=font_size)  # no date / wavelength
    #plt.title(r'%s: %i $\AA$  [Histogram - %s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
    plt.xlabel('%s' % cbar_labels[i], fontsize=font_size, labelpad=10)
    plt.ylabel('Bin Count', fontsize=font_size, labelpad=10)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim(h_min, h_max)
    y, x, _ = plt.hist(flat_param, bins=200, range=(h_min, h_max))
    plt.ylim(0, y.max()*1.1)
    #plt.hist(flatten_slopes, bins='auto')  # try this (actually think we want constant bins throughout wavelengths)
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_Histogram_%s.jpeg' % (directory, date, wavelength, date, wavelength, names[i]))
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_Histogram_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf')
    #"""
