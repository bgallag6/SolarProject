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

#HEATMAPS = np.load('C:/Users/Brendan/Desktop/20130626_171_test/171/param_6_double.npy')
#HEATMAPS2 = np.load('C:/Users/Brendan/Desktop/20130626_171_test/171/param_6_single.npy')
#HEATMAPS3 = np.load('C:/Users/Brendan/Desktop/20130626_171_test/171/param_4_double.npy')
#HEATMAPS4 = np.load('C:/Users/Brendan/Desktop/20130626_171_test/171/param_4_single.npy')
#VISUAL = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_visual.npy')

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20130626'
#h1 = np.load('%s/DATA/Output/%s/171/param.npy' % (directory, date))
#h2 = np.load('%s/DATA/Output/%s/193/param.npy' % (directory, date))
#h3 = np.load('%s/DATA/Output/%s/211/param.npy' % (directory, date))
#h4 = np.load('%s/DATA/Output/%s/304/param.npy' % (directory, date))
h5 = np.load('%s/DATA/Output/%s/1600/param.npy' % (directory, date))
h6 = np.load('%s/DATA/Output/%s/1700/param.npy' % (directory, date))

#r = ss.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20130626', wavelength=193, path_name='C:/Users/Brendan/Desktop/193_midfile_1600x1600_slope6/')
path_name='C:/Users/Brendan/Desktop/histCompare'

# create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
#titles = ['Slope Coefficient', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location [sec]', 'Gaussian Width', '$\chi^2$']
titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'(b) Power Law Index $n$', r'Power Law Tail - $C$', r'Gaussian Amplitude [flux] - $\alpha$', r'Gauss. Loc. $\beta$ [min]', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', 'p-Value']
#names = ['PL_A', 'Slopes', 'PL_C', 'Gauss_Amp', 'Gauss_Loc', 'Gauss_Wid', 'Chi2']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']
#cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '$\chi^2$']
#cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [seconds]', 'Width', '$\chi^2$']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [min]', 'Width', 'F-Statistic', 'Amplitude Scaled', 'P-Value']

#vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore  (or option to set ranges for specific wavelengths?)
#vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore

#M2_low = [0., 0.3, 0., 0.00001, 1./np.exp(-4.6), 0.05]
#M2_low = [-0.0000001, 0.85, -0.0017, 0.00001, (1./np.exp(-4.6))/60., 0.05, 0., 0., 0.]
#M2_high = [0.000035, 1.8, 0.0001, 0.03, (1./np.exp(-6.5))/60., 0.8, 500., 7., 1.]
M2_low = [-0.0000001, 0.85, -0.0017, 0.00001, 1./np.exp(-4.6)/60., 0.05, 0., 0., 0.]
M2_high = [0.000035, 1.8, 0.0001, 0.03, (1./np.exp(-6.5))/60., 0.8, 500., 7., 1.]

wavelength = [1600,1700]
year = date[0:4]
month = date[4:6]
day = date[6:8]
date_title = '%s-%s-%s' % (year,month,day)

# generate p-value heatmap + masked Gaussian component heatmaps
df1, df2 = 3, 6  # degrees of freedom for model M1, M2
p_val5 = ff.sf(h5[6], df1, df2)
p_val6 = ff.sf(h6[6], df1, df2)

mask_thresh = 0.000005  # significance threshold - masked above this value
   
# mask the Gaussian component arrays with NaNs if above threshold
loc_mask5 = np.copy(h5[4])
loc_mask5[p_val5 > mask_thresh] = np.NaN

loc_mask6 = np.copy(h6[4])
loc_mask6[p_val6 > mask_thresh] = np.NaN


heatmap = [h5,h6]
#heatmap = [h1,h2,h3,h4,h5]
#wavelengths = [171,171]
#wavelengths = [193,193]
#wavelengths = [211,211]
#wavelengths = [193 for i in range(len(heatmap))]
#wavelengths = [17100,1930,211,30400,160000]
wavelengths = [1600,1700]

h_map = h5
h_map2 = h6
#h_map = HEATMAPS
#h_map2 = HEATMAPS2
#h_map3 = HEATMAPS3
#h_map4 = HEATMAPS4
#h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]  # trim last row and column from array (originally needed since went one past)
#h_map2 = h_map2[:,0:h_map2.shape[1]-1,0:h_map2.shape[2]-1]  # trim last row and column from array (originally needed since went one past)

#h_map[4] = loc_mask5
#h_map2[4] = loc_mask6


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


#for i in range(0,len(titles)-1):
for i in range(4,5):
    
    
    
    if i == 6:
        NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        #cmap = 'jet'
        cmap = cm.get_cmap('jet', 10)
    elif i == 4:
        h_map[i] = (1./(np.exp(h_map[i])))/60.
        h_map2[i] = (1./(np.exp(h_map2[i])))/60.
        #h_map3[i] = 1./(np.exp(h_map3[i]))
        #h_map4[i] = 1./(np.exp(h_map4[i]))
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
    #fig = plt.figure(figsize=(13,9))
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
    plt.title(r'1600 $\AA$: %s' % (titles[i]), y = 1.01, fontsize=25)  # no date / wavelength
    
    im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=M2_low[i], vmax=M2_high[i])
    #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
    #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    c_ticks = np.zeros((11))
    h_step = (M2_high[i] - M2_low[i])/10
    for h in range(11):
        #c_ticks[h] = h1 + h_step*h
        c_ticks[h] = M2_low[i] + h_step*h
    if i == 0:
        cbar = plt.colorbar(im,cax=cax, format='%0.2e')
    else:
        cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=5) 
    cbar.set_ticks(c_ticks)
    #plt.tight_layout()
    #plt.savefig('%s/%s_%i_heatmap_%s.jpeg' % (path_name, date, wavelength, names[i]))
    plt.savefig('%s/%s_1600_%s.pdf' % (path_name, date, names[i]), format='pdf')
    #"""
    #"""
    
    fig = plt.figure(figsize=(fig_width,fig_height))
    #fig = plt.figure(figsize=(12,10))
    ax = plt.gca()  # get current axis -- to set colorbar 
    #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
    plt.title(r'1700 $\AA$: %s' % (titles[i]), y = 1.01, fontsize=25)  # no date / wavelength
    im = ax.imshow(np.flipud(h_map2[i]), cmap = cmap, vmin=M2_low[i], vmax=M2_high[i])
    #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
    #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    c_ticks = np.zeros((11))
    h_step = (M2_high[i] - M2_low[i])/10
    for h in range(11):
        #c_ticks[h] = h1 + h_step*h
        c_ticks[h] = M2_low[i] + h_step*h
    if i == 0:
        cbar = plt.colorbar(im,cax=cax, format='%0.2e')
    else:
        cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=5) 
    cbar.set_ticks(c_ticks)
    plt.savefig('%s/%s_1700_%s.pdf' % (path_name, date, names[i]), format='pdf')
    
    
    flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
    flat_param2 = np.reshape(h_map2[i], (h_map2[i].shape[0]*h_map2[i].shape[1]))
    #flat_param3 = np.reshape(h_map3[i], (h_map3[i].shape[0]*h_map3[i].shape[1]))
    #flat_param4 = np.reshape(h_map4[i], (h_map4[i].shape[0]*h_map4[i].shape[1]))
         
    # calculate some statistics
    mean = np.mean(flat_param)
    mean2 = np.mean(flat_param2)
    sigma = np.std(flat_param)
    sigma2 = np.std(flat_param2)



    #fig = plt.figure(figsize=(fig_width,fig_height))
    fig = plt.figure(figsize=(13,13))
    #plt.title(r'1600 $\AA$ vs 1700 $\AA$: %s' % titles[i], y = 1.01, fontsize=25)
    plt.xlabel('%s' % cbar_labels[i], fontsize=20, labelpad=10)
    plt.ylabel('Bin Count', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    plt.xlim(M2_low[i],M2_high[i])
    y, x, _ = plt.hist(flat_param, bins=200, color='black', range=(M2_low[i],M2_high[i]), label = r'1600 $\AA$')
      
    n=y[1:-2] 
    bins=x[1:-2]
    elem = np.argmax(n)
    bin_max = bins[elem]
    
    plt.vlines(bin_max, 0, 375000, color='black', linestyle='solid', linewidth=2.) 
    plt.vlines(bin_max, 0, 375000, color='white', linestyle='dashed', linewidth=1.5)       
    #plt.vlines(mean, 0, y.max()*1.1, color='red', linestyle='solid', linewidth=1.5)
    #plt.vlines(0, 0, y.max()*1.1, color='white', linestyle='dashed', linewidth=1.5) 
    
    
    y2, x2, _ = plt.hist(flat_param2, bins=200, color='red', alpha=0.5, range=(M2_low[i],M2_high[i]), label = r'1700 $\AA$')    
    
    n2=y2[1:-2]
    bins2=x2[1:-2]
    elem2 = np.argmax(n2)
    bin_max2 = bins2[elem2]
    
    plt.vlines(bin_max2, 0, 375000, color='red', linestyle='solid', linewidth=2.)
    plt.vlines(bin_max2, 0, 375000, color='white', linestyle='dashed', linewidth=1.5)
    #plt.vlines(mean2, 0, y.max()*1.1, color='red', linestyle='solid', linewidth=1.5)
    #plt.vlines(0, 0, y.max()*1.1, color='white', linestyle='dashed', linewidth=1.5) 
    
    #y, x, _ = plt.hist(flat_param3, bins=200, color='blue', alpha=0.5, range=(M2_low[i],M2_high[i]), label = '4 - Double')
    #y, x, _ = plt.hist(flat_param4, bins=200, color='green', alpha=0.5, range=(M2_low[i],M2_high[i]), label = '4 - Single')
    #plt.ylim(0, y.max()*1.1)
    plt.ylim(0, 375000)
    plt.legend(loc = 'upper right', fontsize=23)
    plt.title('%s \n' % titles[i] + r'1600 $\AA$: mode = %0.3f | mean = %0.3f | Std. Dev. = %0.3f' % (bin_max, mean, sigma) + '\n' + '1700 $\AA$: mode = %0.3f | mean = %0.3f | Std. Dev. = %0.3f' % (bin_max2, mean2, sigma2), y = 1.01, fontsize=21)
    #plt.hist(flatten_slopes, bins='auto')  # try this (actually think we want constant bins throughout wavelengths)
    #plt.savefig('%s/%s_%i_Histogram_%s.jpeg' % (path_name, date, wavelength, names[i]))
    plt.savefig('%s/%s_1600_1700_Histogram_%s.pdf' % (path_name, date, names[i]), format='pdf')
    #"""