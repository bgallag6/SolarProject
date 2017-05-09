# -*- coding: utf-8 -*-
"""
Created on Sat May 06 16:37:39 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm
from scipy import stats
import sunpy
from sunpy.map import Map

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20140818'
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

visual = visual[:,1:-1,1:-1]
#h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]  # trim last row and column from array (originally needed since went one past)

"""
# trim x/y dimensions equally so that resulting region is 1600x1600    
trim_y = (h_map.shape[1]-1600)/2
trim_x = (h_map.shape[2]-1600)/2
h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
x_ind = [-800,-600,-400,-200,0,200,400,600,800]
y_ind = [800,600,400,200,0,-200,-400,-600,-800]    
"""

#h_map = h_map[:, 0:h_map.shape[1]-50, 0:500]  # for 20130626 blobs      
#x_ticks = [0,100,200,300,400,500]
#y_ticks = [0,100,200,300,400]   
#y_ticks = [0,100,200,300] # 20120923

#x_ticks = [0,100,200,300,400,500]
#y_ticks = [0,100,200,300]  


# generate p-value heatmap + masked Gaussian component heatmaps
df1, df2 = 3, 6  # degrees of freedom for model M1, M2
p_val = ff.sf(h_map[6], df1, df2)

p_mask = np.copy(p_val)

mask_thresh = 0.005  # significance threshold - masked above this value
   
p_mask = np.copy(p_val)
p_maskI = np.copy(p_val)

v_mask = np.copy(visual[0])
#v_mask = v_mask[1:-1,1:-1]

# mask the Gaussian component arrays with NaNs if above threshold 
p_mask[p_val > mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
p_maskI[p_val < mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
v_mask[p_val < mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
v_mask[p_val < mask_thresh] = 1.  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN

# determine percentage of region masked 
count = np.count_nonzero(np.isnan(p_mask))   
total_pix = p_val.shape[0]*p_val.shape[1]
mask_percent = ((np.float(count))/total_pix)*100
                
plots = [p_mask, v_mask]  # make array of masked plots to iterate over

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



# generate visual images
titles_vis = ['Average', 'Middle-File']
names_vis = ['average', 'mid']

vis = visual
#trim_yv = (vis.shape[1]-1600)/2
#trim_xv = (vis.shape[2]-1600)/2
#vis = vis[:, trim_yv:vis.shape[1]-trim_yv, trim_xv:vis.shape[2]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

#vis = vis[:, 0:vis.shape[1]-50, 0:500]  # for 20130626 blobs     

for i in range(1):
    
    v_min = np.percentile(vis[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    v_max = np.percentile(vis[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)     
    
    delta = 1.
    x = np.arange(0., v_mask.shape[1], delta)
    y = np.arange(0., v_mask.shape[0], delta)
    X, Y = np.meshgrid(x, y)
    #Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    #Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    # difference of Gaussians
    #Z = 10.0 * (Z2 - Z1)
    Z = np.flipud(v_mask)
    
    #fig = plt.figure(figsize=(12,9))
    fig = plt.figure(figsize=(fig_width,fig_height))
    
    ax = plt.gca()
    #plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
    plt.title('Visual: %s' % (titles_vis[i]), y = 1.02, fontsize=font_size)  # no date / wavelength
    #im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
    #im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
    #cmap = cm.get_cmap('sdoaia%i' % wavelength, 10)    
    #im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
    im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
    CS = plt.contour(X, Y, Z, levels=[2.], linewidths=2, colors='white', linestyles='solid')
    CS = plt.contour(X, Y, Z, levels=[2.], linewidths=2, colors='black', linestyles='dashed')
    #im = ax.imshow(np.flipud(v_mask), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
    #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
    #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
    #plt.xticks(x_ticks,fontsize=font_size)
    #plt.yticks(y_ticks,fontsize=font_size)
    #plt.xticks(fontsize=font_size)
    #plt.yticks(fontsize=font_size)
    #plt.xticks(x_ticks,x_ind,fontsize=font_size)
    #plt.yticks(y_ticks,y_ind,fontsize=font_size)
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('Intensity', size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    #plt.tight_layout()
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
    #plt.savefig('C:/Users/Brendan/Desktop/pvalue_mask_contour_overlay.pdf', format='pdf')
    
v_min = np.percentile(vis[0],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
v_max = np.percentile(vis[0],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)     
delta = 1.
x = np.arange(0., v_mask.shape[1], delta)
y = np.arange(0., v_mask.shape[0], delta)
X, Y = np.meshgrid(x, y)
#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
#Z = 10.0 * (Z2 - Z1)
Z = np.flipud(v_mask)

#fig = plt.figure(figsize=(12,9))
fig = plt.figure(figsize=(fig_width,fig_height))

ax = plt.gca()
#plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
plt.title('Visual: %s' % (titles_vis[0]), y = 1.02, fontsize=font_size)  # no date / wavelength
#im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
#im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#cmap = cm.get_cmap('sdoaia%i' % wavelength, 10)    
#im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
im = ax.imshow(np.flipud(v_mask), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
CS = plt.contour(X, Y, Z, levels=[2.], linewidths=2, colors='white', linestyles='solid')
CS = plt.contour(X, Y, Z, levels=[2.], linewidths=2, colors='black', linestyles='dashed')
#plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
#plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
#plt.xticks(x_ticks,fontsize=font_size)
#plt.yticks(y_ticks,fontsize=font_size)
#plt.xticks(fontsize=font_size)
#plt.yticks(fontsize=font_size)
#plt.xticks(x_ticks,x_ind,fontsize=font_size)
#plt.yticks(y_ticks,y_ind,fontsize=font_size)
ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
#cbar.set_label('Intensity', size=20, labelpad=10)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
#plt.tight_layout()
#plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
#plt.savefig('C:/Users/Brendan/Desktop/pvalue_mask_inverted_create_contour.pdf', format='pdf')



h_min = 0.0
h_max = 0.005
#fig = plt.figure(figsize=(12,9))
fig = plt.figure(figsize=(fig_width,fig_height))

ax = plt.gca()
#plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
plt.title('p-Value Mask', y = 1.02, fontsize=font_size)  # no date / wavelength
#im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
#im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#cmap = cm.get_cmap('sdoaia%i' % wavelength, 10)    
#im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
im = ax.imshow(np.flipud(p_mask), vmin = h_min, vmax = h_max)
ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
#cbar.set_label('Intensity', size=20, labelpad=10)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
#plt.tight_layout()
#plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
#plt.savefig('C:/Users/Brendan/Desktop/pvalue_mask_inverted.pdf', format='pdf')

h_min = 0.0
h_max = 0.005
#fig = plt.figure(figsize=(12,9))
fig = plt.figure(figsize=(fig_width,fig_height))

ax = plt.gca()
#plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
plt.title('p-Value Mask', y = 1.02, fontsize=font_size)  # no date / wavelength
#im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
#im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#cmap = cm.get_cmap('sdoaia%i' % wavelength, 10)    
#im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
im = ax.imshow(np.flipud(p_maskI), vmin = h_min, vmax = h_max)
ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
#cbar.set_label('Intensity', size=20, labelpad=10)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
#plt.tight_layout()
#plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
#plt.savefig('C:/Users/Brendan/Desktop/pvalue_mask_normal.pdf', format='pdf')