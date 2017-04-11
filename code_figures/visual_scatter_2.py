# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 03:25:21 2017

@author: Brendan
"""
        
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
import numpy as np
import astropy.units as u
import matplotlib.patches as patches

v171 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/visual.npy')
p171 = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/171/param.npy')
#v171 = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_171_-500_500i_-500_600j_visual.npy')
#p171 = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_171_-500_500i_-500_600j_param_slope6_arthm.npy')
#v1600 = np.load('C:/Users/Brendan/Desktop/1600/visual_1600.npy')


    
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


#m2 = [743, 708, 525, 757, 765, 722, 867]
#width = np.array([v171[0].shape[1] for i in range(len(m2))])
#m2 = width - m2
#l2 = [322, 352, 551, 319, 325, 1441, 864]
#height = np.array([v171[0].shape[0] for i in range(len(l2))])
#l2 = height - l2

#[743]
#[322]

m2 = [722, 188, 726, 872]
l2 = [v171[0].shape[0]-1427, v171[0].shape[0]-523, v171[0].shape[0]-328, v171[0].shape[0]-875]

m3 = [100]
l3 = [1600-900]

wavelength = 171

date_title = '2013/06/26'
h_map = v171


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


  
# generate visual images
titles_vis = ['Average', 'Middle-File']
names_vis = ['average', 'mid']

vis = v171
trim_yv = (vis.shape[1]-1600)/2
trim_xv = (vis.shape[2]-1600)/2
#vis = vis[:, trim_yv:vis.shape[1]-trim_yv, trim_xv:vis.shape[2]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)  
#vis1600 = v1600[:, 19:v1600.shape[1]-19, 12:v1600.shape[2]-11]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)  

#vis1600_trim = vis1600[:, :, 1100:1400] 

for i in range(1):
    
    v_min = np.percentile(vis[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    v_max = np.percentile(vis[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)     
    
    #v_minB = np.percentile(vis1600[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    #v_maxB = np.percentile(vis1600[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)    
    
    """
    #fig = plt.figure(figsize=(12,9))
    fig = plt.figure(figsize=(fig_width,fig_height))
    
    vis[i,800:1000,1000:1300] = np.NaN  
    #vis = np.ma.masked_where(vis == 0, vis)
    
    ax = plt.gca()
    #plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
    plt.title(r'171$\AA$  Visual: %s w/ Points' % (titles_vis[i]), y = 1.01, fontsize=23)  # no date / wavelength
    #im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
    im = ax.imshow(np.flipud(vis304[i]), cmap='sdoaia%i' % 304, vmin = 0, vmax = 500)
    im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max) 
    
    plt.scatter(m2, l2, s=75, c='white')
    #plt.scatter(m2, l2, s=75, c='red', marker='*')
    plt.xlim(0, vis[0].shape[1])
    plt.ylim(vis[0].shape[0], 0)
    plt.xlabel('X-Position [Pixels]', fontsize=23, labelpad=10)
    plt.ylabel('Y-Position [Pixels]', fontsize=23, labelpad=10)
    plt.xticks([0,200,400,600,800,1000,1200,1400,1600],fontsize=23)
    plt.yticks([0,200,400,600,800,1000,1200,1400,1600],fontsize=23)
    rect = patches.Rectangle((610,1600-1460), 70, 90, color='white', fill=True)
    ax.add_patch(rect)
    ax.text(620,1600-1475+90, 'A', fontsize=25)
    rect = patches.Rectangle((580,1600-580), 70, 90, color='white', fill=True)
    ax.add_patch(rect)
    ax.text(590,1600-595+90, 'B', fontsize=25)   
    rect = patches.Rectangle((815,1600-330), 70, 90, color='white', fill=True)
    ax.add_patch(rect)
    ax.text(825,1600-345+90, 'C', fontsize=25)
    rect = patches.Rectangle((1190,570), 70, 90, color='white', fill=True)
    ax.add_patch(rect)
    ax.text(1200,645, 'D', fontsize=25)
    
    
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    
    #cbar.set_label('Intensity', size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=23, pad=5) 
    #plt.tight_layout()
    #plt.savefig('%s/%s_%i_visual_%s.jpeg' % (path_name, date, wavelength, names_vis[i]))
    #plt.savefig('C:/Users/Brendan/Desktop/171_visual_%s_points_overlay.pdf' % names_vis[i], format='pdf')
    """
    
    x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
    y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
    x_ind = [-800,-600,-400,-200,0,200,400,600,800]
    y_ind = [800,600,400,200,0,-200,-400,-600,-800]
   

    fig = plt.figure(figsize=(15,15))
    #ax1 = plt.subplot2grid((1,6),(0, 0), colspan=5, rowspan=1)  
    ax1 = plt.gca()
    plt.title(r'171$\AA$  Visual: %s w/ Points' % (titles_vis[i]), y = 1.02, fontsize=23)  # no date / wavelength
    im = ax1.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
    plt.scatter(m2, l2, s=100, c='red')
    plt.scatter(m2, l2, s=50, c='white')
    #plt.scatter(m2, l2, s=75, c='red', marker='*')
    plt.xlim(0, vis[0].shape[1])
    plt.ylim(vis[0].shape[0], 0)
    #plt.xlabel('X-Position [Pixels]', fontsize=23, labelpad=10)
    #plt.ylabel('Y-Position [Pixels]', fontsize=23, labelpad=10)
    #plt.xticks([0,200,400,600,800,1000,1200,1400,1600],fontsize=23)
    #plt.yticks([0,200,400,600,800,1000,1200,1400,1600],fontsize=23)
    plt.xticks(x_ticks,x_ind,fontsize=23)
    plt.yticks(y_ticks,y_ind,fontsize=23)
    ax1.tick_params(axis='both', which='major', pad=10)
    rect = patches.Rectangle((610,140), 70, 90, color='white', fill=True)
    ax1.add_patch(rect)
    ax1.text(623.2,209, 'A', fontsize=27)
    rect = patches.Rectangle((220,1160), 70, 90, color='white', fill=True)
    ax1.add_patch(rect)
    ax1.text(231,1229.5, 'C', fontsize=27)   
    rect = patches.Rectangle((770,1230), 70, 90, color='white', fill=True)
    ax1.add_patch(rect)
    ax1.text(782.2,1298.7, 'D', fontsize=27)
    rect = patches.Rectangle((905,810), 70, 90, color='white', fill=True)
    ax1.add_patch(rect)
    ax1.text(918.3,877.6, 'B', fontsize=27)
    #rect6 = patches.Rectangle((1100,1), 300, 1600, color='white', fill=False, linewidth=3)
    #ax1.add_patch(rect6)
    #ax1.text(1180, 110, ' Inset\nRegion', fontsize=21, color='white')
    rect = patches.Rectangle((906,811), 68, 88, color='red', fill=False, linewidth=2)
    ax1.add_patch(rect)
    rect = patches.Rectangle((221,1161), 68, 88, color='red', fill=False, linewidth=2)
    ax1.add_patch(rect)
    rect = patches.Rectangle((771,1231), 68, 88, color='red', fill=False, linewidth=2)
    ax1.add_patch(rect)
    rect = patches.Rectangle((611,141), 68, 88, color='red', fill=False, linewidth=2)
    ax1.add_patch(rect)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('Intensity', size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=23, pad=5)

    """
    ax2 = plt.subplot2grid((1,6),(0, 5), colspan=1, rowspan=1)  
    plt.title(r'1600$\AA$', y = 1.01, fontsize=23)  # no date / wavelength
    im = ax2.imshow(np.flipud(vis1600_trim[i]), cmap='sdoaia1600', vmin = v_minB, vmax = v_maxB)
    plt.scatter(m3, l3, s=75, c='white')
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.xlim(0, vis1600_trim[0].shape[1])
    plt.ylim(vis1600_trim[0].shape[0], 0)
    plt.xlabel('X-Position [Pixels]', fontsize=23, labelpad=10, visible=False)
    rect8 = patches.Rectangle((130,570), 70, 90, color='white', fill=True)
    ax2.add_patch(rect8)
    ax2.text(143,640, 'D', fontsize=25)
    plt.subplots_adjust(wspace=0)
    ax2.text(90, 60, 'Inset', fontsize=21, color='white')
    rect9 = patches.Rectangle((1,1), 298, 1598, color='white', fill=False, linewidth=3)
    ax2.add_patch(rect9)
    """
    plt.savefig('C:/Users/Brendan/Desktop/171_final_points/171_visual_%s_points_overlay_revised7.jpeg' % names_vis[i])
    #plt.savefig('C:/Users/Brendan/Desktop/171_final_points/171_visual_%s_points_overlay_revised6.pdf' % names_vis[i], format='pdf')