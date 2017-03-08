# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 17:01:57 2017

@author: Brendan
"""


import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import f as ff

#HEATMAPS = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20130530_1600_2300_2600_2200_3000_float_numpy.npy')
HEATMAPS = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_171_-500_500i_-500_600j_param_slope6_arthm.npy')


titles = ['Power Law Slope Coefficient', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width', '($/chi^2$)']
names = ['PL_A', 'Slopes', 'PL_C', 'Gauss_Amp', 'Gauss_Loc', 'Gauss_Wid', 'Chi2']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '($/chi^2$)']
#vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore
#vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore
#wavelength = wavelength
#year = date[0:4]
#month = date[4:6]
#day = date[6:8]
wavelength = 171
year = '2013'
month = '06'
day = '26'
date_title = '%s-%s-%s' % (year,month,day)

h_map = HEATMAPS


# generate p-value heatmap
df1, df2 = 3, 6
p_val = ff.sf(h_map[6], df1, df2)

p_mask = np.copy(p_val)

mask_arr = [0.005,0.0005,0.00005,0.000001]
for j in range(len(mask_arr)):
    mask_thresh = mask_arr[j]
       
    p_mask = np.copy(p_val)
    amp_mask = np.copy(h_map[3])
    loc_mask = np.copy(h_map[4])
    wid_mask = np.copy(h_map[5])
    
    h_min_amp = np.percentile(h_map[3],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max_amp = np.percentile(h_map[3],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    
    count = 0
    
    for i in range(p_val.shape[0]):
            for j in range(p_val.shape[1]):
                if p_val[i][j] > mask_thresh:
                    count += 1
                    p_mask[i][j] = np.NaN
                    amp_mask[i][j] = np.NaN
                    loc_mask[i][j] = np.NaN
                    wid_mask[i][j] = np.NaN
    
    flat_hmap_amp = np.reshape(amp_mask,(h_map.shape[1]*h_map.shape[2]))
    flat_hmap_loc = np.reshape(loc_mask,(h_map.shape[1]*h_map.shape[2]))
    flat_hmap_wid = np.reshape(wid_mask,(h_map.shape[1]*h_map.shape[2]))
    flat_hmap_index = np.reshape(h_map[1],(h_map.shape[1]*h_map.shape[2]))
    #flat_hmap = np.zeros((h_map.shape[0],h_map.shape[1]*h_map.shape[2]))

    fig = plt.figure(figsize=(15,9))
    xmin = np.percentile(flat_hmap_amp,1)
    xmax = np.percentile(flat_hmap_amp,99)
    ymin = np.percentile(flat_hmap_index,1)
    ymax = np.percentile(flat_hmap_index,99)
    x_margin = (xmax-xmin)*0.05
    y_margin = (ymax-ymin)*0.05
    
    #ax = fig.add_subplot(111,projection='3d')
    ax = fig.add_subplot(111)
    ax.set_title('Mask Threshold = %f' % mask_thresh)
    #ax.scatter(fl_plC, fl_slopes, fl_plA,marker='.')
    ax.set_xlim(-0.005, 0.04)
    ax.set_ylim(0.75, 3.25)
    ax.set_xlabel('Gaussian Amplitude')
    ax.set_ylabel('Power Law Index')
    ax.scatter(flat_hmap_amp,flat_hmap_index,marker='.')
    plt.savefig('C:/Users/Brendan/Desktop/171_scatter_mask_%i_amp.pdf' % mask_thresh, format='pdf')
    """
    fig = plt.figure(figsize=(15,9))
    xmin = np.percentile(flat_hmap_loc,1)
    xmax = np.percentile(flat_hmap_loc,99)
    ymin = np.percentile(flat_hmap_index,1)
    ymax = np.percentile(flat_hmap_index,99)
    x_margin = (xmax-xmin)*0.05
    y_margin = (ymax-ymin)*0.05
    
    #ax = fig.add_subplot(111,projection='3d')
    ax = fig.add_subplot(111)
    #ax.scatter(fl_plC, fl_slopes, fl_plA,marker='.')
    ax.set_xlim(-6.5,-4.6)
    ax.set_ylim(0.75, 3.25)
    ax.set_xlabel('Gaussian Location')
    ax.set_ylabel('Power Law Index')
    ax.scatter(flat_hmap_loc,flat_hmap_index,marker='.')
    #plt.savefig('C:/Users/Brendan/Desktop/171_scatter_mask_%i_loc.pdf' % mask_thresh, format='pdf')
    
    fig = plt.figure(figsize=(15,9))
    xmin = np.percentile(flat_hmap_wid,1)
    xmax = np.percentile(flat_hmap_wid,99)
    ymin = np.percentile(flat_hmap_index,1)
    ymax = np.percentile(flat_hmap_index,99)
    x_margin = (xmax-xmin)*0.05
    y_margin = (ymax-ymin)*0.05
    
    #ax = fig.add_subplot(111,projection='3d')
    ax = fig.add_subplot(111)
    #ax.scatter(fl_plC, fl_slopes, fl_plA,marker='.')
    ax.set_xlim(0.05,0.8)
    ax.set_ylim(0.75, 3.25)
    ax.set_xlabel('Gaussian Width')
    ax.set_ylabel('Power Law Index')
    ax.scatter(flat_hmap_wid,flat_hmap_index,marker='.')
    #plt.savefig('C:/Users/Brendan/Desktop/171_scatter_mask_%i_wid.pdf' % mask_thresh, format='pdf')
    """
"""
for i in range(h_map.shape[0]):
    if i == 6:
        NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        flat_hmap[i] = np.reshape(NaN_replace, (h_map.shape[1]*h_map.shape[2]))
    else:
        flat_hmap[i] = np.reshape(h_map[i], (h_map.shape[1]*h_map.shape[2]))


#for j in range(h_map.shape[0]):
for j in range(1,2):
    #for i in range(h_map.shape[0]): 
    for i in range(3,4): 
        if i > j:
            xmin = np.percentile(flat_hmap[j],1)
            xmax = np.percentile(flat_hmap[j],99)
            ymin = np.percentile(flat_hmap[i],1)
            ymax = np.percentile(flat_hmap[i],99)
            x_margin = (xmax-xmin)*0.05
            y_margin = (ymax-ymin)*0.05
            fig = plt.figure(figsize=(15,9))
            #ax = fig.add_subplot(111,projection='3d')
            ax = fig.add_subplot(111)
            #ax.scatter(fl_plC, fl_slopes, fl_plA,marker='.')
            ax.set_xlim(xmin-x_margin, xmax+x_margin)
            ax.set_ylim(ymin-y_margin, ymax+y_margin)
            ax.set_xlabel('%s' % titles[j])
            ax.set_ylabel('%s' % titles[i])
            ax.scatter(flat_hmap[j],flat_hmap[i],marker='.')
            #path_name='C:/Users/Brendan/Desktop/PHYS 326/test_temp'
            #date = '20130530'
            #plt.savefig('%s/%s_%i_scatter_%s_vs_%s.jpeg' % (path_name, date, wavelength, names[i], names[j]))
"""  
    
"""
for i in range(0,len(titles)):    
    fig = plt.figure(figsize=(15,9))
    ax = plt.gca()
    plt.title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[i]), y = 1.01, fontsize=25)
    
    if i == 6:
        NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    else:
        h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    
    im = ax.imshow(h_map[i], vmin=h_min, vmax=h_max)
    plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
    plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=5) 
    plt.tight_layout()
    #plt.savefig('%s/%s_%i_heatmap_%s.jpeg' % (path_name, date, wavelength, names[i]))
    plt.savefig('%s/%s_%i_heatmap_%s.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
"""   
    
