# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 14:17:51 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm

HEATMAPS = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_param.npy')
HEATMAPS2 = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_param_slope6.npy')
VISUAL = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_visual.npy')

#r = ss.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20130626', wavelength=193, path_name='C:/Users/Brendan/Desktop/193_midfile_1600x1600_slope6/')
date = '20130626'
wavelength=193
path_name='C:/Users/Brendan/Desktop/193_midfile_1600x1600_comparison/'

# create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
#titles = ['Slope Coefficient', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location [sec]', 'Gaussian Width', '$\chi^2$']
titles = [r'Power Law Slope-Coefficient -- [$A$]', r'Power Law Index -- [$n$]', r'Power Law Tail -- [$C$]', r'Gaussian Amplitude -- [$\alpha$]', r'Gaussian Location [Seconds] -- [$\beta$]', r'Gaussian Width -- [$\sigma$]', 'F-Statistic', r'Gaussian Amplitude Scaled -- [$\alpha$]', 'P-Value']
#names = ['PL_A', 'Slopes', 'PL_C', 'Gauss_Amp', 'Gauss_Loc', 'Gauss_Wid', 'Chi2']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']
#cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '$\chi^2$']
#cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [seconds]', 'Width', '$\chi^2$']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [seconds]', 'Width', 'F-Statistic', 'Amplitude Scaled', 'P-Value']

#vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore  (or option to set ranges for specific wavelengths?)
#vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore

M2_low = [0., 0.3, 0., 0.00001, 1./np.exp(-4.6), 0.05]
M2_high = [0.000002, 3.0, 0.003, 0.01, 1./np.exp(-6.5), 0.8]

wavelength = wavelength
year = date[0:4]
month = date[4:6]
day = date[6:8]
date_title = '%s-%s-%s' % (year,month,day)

h_map = HEATMAPS
h_map2 = HEATMAPS2
h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]  # trim last row and column from array (originally needed since went one past)
h_map2 = h_map2[:,0:h_map2.shape[1]-1,0:h_map2.shape[2]-1]  # trim last row and column from array (originally needed since went one past)

if h_map.shape[2] > h_map.shape[1]:
    aspect_ratio = float(h_map.shape[2]) / float(h_map.shape[1])
    fig_height = 10
    fig_width = 10*aspect_ratio
    
else:
    aspect_ratio = float(h_map.shape[1]) / float(h_map.shape[2])
    #print aspect_ratio
    #fig_width = 10
    fig_width = 10+1  # works better for 20130626 (with no x/y labels)
    #fig_height = 10*aspect_ratio
    fig_height = 10*aspect_ratio  # works better for 20130626


for i in range(0,len(titles)-1):
#for i in range(0,1):
    
    #fig = plt.figure(figsize=(13,9))
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
    plt.title('%s - Index Bounds [0.3, 4.0]' % (titles[i]), y = 1.01, fontsize=25)  # no date / wavelength
    
    if i == 6:
        NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        #cmap = 'jet'
        cmap = cm.get_cmap('jet', 10)
    elif i == 4:
        h_map[i] = 1./(np.exp(h_map[i]))
        h_map2[i] = 1./(np.exp(h_map2[i]))
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
    #"""
    im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=M2_low[i], vmax=M2_high[i])
    #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
    #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    if i == 0:
        cbar = plt.colorbar(im,cax=cax, format='%0.2e')
    else:
        cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=5) 
    #plt.tight_layout()
    #plt.savefig('%s/%s_%i_heatmap_%s.jpeg' % (path_name, date, wavelength, names[i]))
    #plt.savefig('%s/%s_%i_%s_4slope.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
    #"""
    #"""
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
    plt.title('%s - Index Bounds [0.3, 6.0]' % (titles[i]), y = 1.01, fontsize=25)  # no date / wavelength
    im = ax.imshow(np.flipud(h_map2[i]), cmap = cmap, vmin=M2_low[i], vmax=M2_high[i])
    #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
    #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    if i == 0:
        cbar = plt.colorbar(im,cax=cax, format='%0.2e')
    else:
        cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=5) 
    plt.savefig('%s/%s_%i_%s_6slope.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
    
    flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
    flat_param2 = np.reshape(h_map2[i], (h_map2[i].shape[0]*h_map2[i].shape[1]))

    fig = plt.figure(figsize=(12,9))
    plt.title('Index Bound Comparison: %s' % titles[i], y = 1.01, fontsize=25)
    plt.xlabel('%s' % cbar_labels[i], fontsize=20, labelpad=10)
    plt.ylabel('Bin Count', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    plt.xlim(M2_low[i],M2_high[i])
    y, x, _ = plt.hist(flat_param, bins=200, color='black', range=(M2_low[i],M2_high[i]), label = '[0.3, 4.0]')
    y, x, _ = plt.hist(flat_param2, bins=200, color='red', alpha=0.5, range=(M2_low[i],M2_high[i]), label = '[0.3, 6.0]')
    plt.ylim(0, y.max()*1.1)
    plt.legend(loc = 'upper left')
    #plt.hist(flatten_slopes, bins='auto')  # try this (actually think we want constant bins throughout wavelengths)
    #plt.savefig('%s/%s_%i_Histogram_%s.jpeg' % (path_name, date, wavelength, names[i]))
    #plt.savefig('%s/%s_%i_Histogram_%s.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
    #"""