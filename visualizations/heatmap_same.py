# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 06:50:03 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm


#h171 = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/171/20130626_171_-500_500i_-500_500j_param.npy')
#h193 = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/193/20130626_193_-500_500i_-500_500j_param.npy')
#h211 = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/211/20130626_211_-500_500i_-500_500j_param.npy')
#h304 = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/304/20130626_304_-500_500i_-500_500j_param.npy')

#h171 = np.load('C:/Users/Brendan/Desktop/project_files/param.npy')
#h171 = np.load('F:/Users/Brendan/Desktop/SolarProject/data_sort/20130626/171/20130626_171_-500_500i_-500_500j_param_newsigma.npy')
#h_old = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_304_-500_500i_-500_500j_param.npy')
#h_new = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_304_-500_500i_-500_600j_param_slope6_arthm.npy')


h1 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param1.npy')
h2 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param2.npy')
h3 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param3.npy')
h4 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param4.npy')
h5 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param5.npy')
h6 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param6.npy')
h7 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param7.npy')
h8 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param8.npy')

    
date = '20130626'
#path_name = 'F:/Users/Brendan/Desktop/SolarProject/data/20130626'
path_name = 'C:/Users/Brendan/Desktop/171_segments'

# create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
#titles = ['Slope Coefficient', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location [sec]', 'Gaussian Width', '$\chi^2$']
titles = [r'Power Law Slope-Coefficient -- [$A$]', r'Power Law Index -- [$n$]', r'Power Law Tail -- [$C$]', r'Gaussian Amplitude -- [$\alpha$]', r'Gaussian Location [Seconds] -- [$\beta$]', r'Gaussian Width -- [$\sigma$]', 'F-Statistic', r'Gaussian Amplitude Scaled -- [$\alpha$]', 'P-Value']
#names = ['PL_A', 'Slopes', 'PL_C', 'Gauss_Amp', 'Gauss_Loc', 'Gauss_Wid', 'Chi2']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']
#cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '$\chi^2$']
#cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [seconds]', 'Width', '$\chi^2$']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [seconds]', 'Width', 'F-Statistic', 'Amplitude Scaled', 'P-Value']



#heatmap = [h171,h193,h211,h304]
#wavelengths = [171,193,211,304]

#heatmap = [h171,h171_new]
heatmap = [h1,h2,h3,h4,h5,h6,h7,h8]
#wavelengths = [171,171]
#wavelengths = [193,193]
#wavelengths = [211,211]
wavelengths = [171,171,171,171,171,171,171,171]

#h_map2 = h_new

year = date[0:4]
month = date[4:6]
day = date[6:8]
date_title = '%s-%s-%s' % (year,month,day)


vmin = np.zeros((7))
vmax = np.zeros((7))

"""
for m in range(6):
    for n in range(4):
        vmin_temp[n] = np.min(heatmap[n][m])
        vmax_temp[n] = np.max(heatmap[n][m])
    vmin[m] = np.min(vmin_temp)
    vmax[m] = np.max(vmax_temp)
"""
#vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore  (or option to set ranges for specific wavelengths?)
#vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore
M2_low = [0., 0.3, 0., 0.00001, 1./np.exp(-4.6), 0.05]
M2_high = [0.000002, 3.0, 0.003, 0.01, 1./np.exp(-6.5), 0.8]

for m in range(7):
        
    for n in range(len(heatmap)):
        if n == 0:
            v_temp = []
        v_temp_flat = np.reshape(heatmap[n][m], (heatmap[n][m].shape[0]*heatmap[n][m].shape[1]))
        v_temp = np.append(v_temp, v_temp_flat)
    if m == 6:
        NaN_replace = np.nan_to_num(v_temp)  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        vmin[m] = np.percentile(NaN_replace,1)
        vmax[m] = np.percentile(NaN_replace,99)
    elif m == 4:
        v_temp_loc = 1./(np.exp(v_temp))
        vmin[m] = np.percentile(v_temp_loc,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        vmax[m] = np.percentile(v_temp_loc,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        #cmap = 'jet_r'  # reverse color-scale for Gaussian Location, because of flipped frequencies to seconds
        cmap = cm.get_cmap('jet_r', 10)
    else:
        vmin[m] = np.percentile(v_temp,3)
        vmax[m] = np.percentile(v_temp,97)



for c in range(len(heatmap)):
    h_map = heatmap[c]
    wavelength = wavelengths[c]

    h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]  # trim last row and column from array (originally needed since went one past)
    
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
    
    
    for i in range(7):
    #for i in range(0,2):
        
        
        if i == 6:
            NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
            h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #cmap = 'jet'
            cmap = cm.get_cmap('jet', 10)
        elif i == 4:
            h_map[i] = 1./(np.exp(h_map[i]))
            #h_map2[i] = 1./(np.exp(h_map2[i]))
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
        """    
        #fig = plt.figure(figsize=(13,9))
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
        plt.title('%s' % (titles[i]), y = 1.01, fontsize=25)  # no date / wavelength
        
        im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=h_min, vmax=h_max)
        #im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=M2_low[i], vmax=M2_high[i])  #other script used this?
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
        plt.savefig('%s/%s_%i_%s_diff_%i.pdf' % (path_name, date, wavelength, names[i], c), format='pdf')
        """
        #fig = plt.figure(figsize=(13,9))
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
        plt.title('%s -- Segment %i' % (titles[i], c), y = 1.01, fontsize=25)  # no date / wavelength
        
        im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=vmin[i], vmax=vmax[i])
        #im = ax.imshow(np.flipud(h_map2[i]), cmap = cmap, vmin=M2_low[i], vmax=M2_high[i])
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
        plt.savefig('%s/%s_%i_%s_same_%i.jpeg' % (path_name, date, wavelength, names[i], c))
        #plt.savefig('%s/%s_%i_%s_same_%i.pdf' % (path_name, date, wavelength, names[i], c), format='pdf')
        """
        if c == 0:
            if i != 6:
                
                
                if i == 4: 
                    flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
                    flat_param2 = 1./np.exp(np.reshape(h_map2[i], (h_map2[i].shape[0]*h_map2[i].shape[1])))
                else:
                    flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
                    flat_param2 = np.reshape(h_map2[i], (h_map2[i].shape[0]*h_map2[i].shape[1]))
            
                fig = plt.figure(figsize=(12,9))
                plt.title('Geometric vs Arithmetic Average: %s' % titles[i], y = 1.01, fontsize=25)
                plt.title(r'%s: %i $\AA$  [Histogram - %s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
                plt.xlabel('%s' % cbar_labels[i], fontsize=20, labelpad=10)
                plt.ylabel('Bin Count', fontsize=20, labelpad=10)
                plt.xticks(fontsize=17)
                plt.yticks(fontsize=17)
                plt.xlim(M2_low[i],M2_high[i])
                y, x, _ = plt.hist(flat_param, bins=200, color='black', range=(M2_low[i],M2_high[i]), label = 'Geometric')
                y, x, _ = plt.hist(flat_param2, bins=200, color='red', alpha=0.5, range=(M2_low[i],M2_high[i]), label = 'Arithmetic')
                plt.ylim(0, y.max()*1.1)
                plt.legend(loc = 'upper left')
                #plt.xlim(h_min, h_max)
                #y, x, _ = plt.hist(flat_param, bins=200, range=(h_min, h_max))
                #plt.hist(flatten_slopes, bins='auto')  # try this (actually think we want constant bins throughout wavelengths)
                #plt.savefig('%s/%s_%i_Histogram_%s.jpeg' % (path_name, date, wavelength, names[i]))
                plt.savefig('%s/%s_%i_Histogram_%s.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
        """
        
    """
    # generate p-value heatmap
    df1, df2 = 3, 6
    p_val = ff.sf(h_map[6], df1, df2)
    
    p_mask = np.copy(p_val)
    
    mask_thresh = 0.005    
       
    p_mask = np.copy(p_val)
    amp_mask = np.copy(h_map[3])
    loc_mask = np.copy(h_map[4])
    wid_mask = np.copy(h_map[5])
    
    h_min_amp = np.percentile(h_map[3],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max_amp = np.percentile(h_map[3],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    
    for i in range(p_val.shape[0]):
            for j in range(p_val.shape[1]):
                if p_val[i][j] > mask_thresh:
                    p_mask[i][j] = np.NaN
                    amp_mask[i][j] = np.NaN
                    loc_mask[i][j] = np.NaN
                    wid_mask[i][j] = np.NaN
                    
    plots = [p_mask, amp_mask, loc_mask, wid_mask]
    names_f = ['P-Value Mask', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width']
    names_m = ['p_mask', 'amp', 'loc', 'wid']
    
   
    for k in range(4):           
    
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        if k == 0:
            plt.title('P-Value < %0.3f' % (mask_thresh), y = 1.01, fontsize=25)
        else:
            plt.title('%s: P-Value < %0.3f' % (names_f[k], mask_thresh), y = 1.01, fontsize=25)
        if k == 2:
            #cmap = 'jet_r'
            cmap = cm.get_cmap('jet_r', 10)
        else:
            #cmap = 'jet'
            cmap = cm.get_cmap('jet', 10)
        if k == 1:
            im = ax.imshow(np.flipud(plots[k]), cmap = cmap, vmin = h_min_amp, vmax = h_max_amp)
        else:
            im = ax.imshow(np.flipud(plots[k]), cmap = cmap)
        #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
        #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=17, pad=5) 
        #plt.tight_layout()
        #plt.savefig('%s/%s_%i_heatmap_%s.jpeg' % (path_name, date, wavelength, names[i]))
        #plt.savefig('%s/%s_%i_%s_mask_%i.pdf' % (path_name, date, wavelength, names_m[k], (1./mask_thresh)), format='pdf')
        
    
    # generate 'rollover frequency' heatmap
    roll_freq = (h_map[2] / h_map[0])**(-1./ h_map[1])
    roll_freq = 1./roll_freq
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
    plt.title(r'Rollover Period [Seconds] -- [$(C/A)^{-\frac{1}{n}}$]', y = 1.01, fontsize=25)  # no date / wavelength
    roll_freq = np.nan_to_num(roll_freq)  # deal with NaN's causing issues
    h_min = np.percentile(roll_freq,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(roll_freq,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    #cmap = 'jet'      
    cmap = cm.get_cmap('jet', 10)    
    
    im = ax.imshow(np.flipud(roll_freq), cmap = cmap, vmin=h_min, vmax=h_max)
    #im = ax.imshow(np.flipud(roll_freq), cmap = cmap, vmin=(1./10**-1.), vmax=(1./10**-3.5))  # should bounds be set at frequency range
    #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
    #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=5) 
    #plt.tight_layout()
    #plt.savefig('%s/%s_%i_heatmap_%s.jpeg' % (path_name, date, wavelength, names[i]))
    #plt.savefig('%s/%s_%i_roll_freq.pdf' % (path_name, date, wavelength), format='pdf')
    """
    
  
