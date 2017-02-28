# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 15:49:45 2017

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


v171 = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/171/20130626_171_-500_500i_-500_500j_visual.npy')
    
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


m2 = [743, 708, 525, 757, 765, 722, 867]
width = np.array([v171[0].shape[1] for i in range(len(m2))])
#m2 = width - m2
l2 = [322, 352, 551, 319, 325, 1441, 864]
height = np.array([v171[0].shape[0] for i in range(len(l2))])
#l2 = height - l2

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

for i in range(2):
    
    v_min = np.percentile(vis[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    v_max = np.percentile(vis[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)     
    
    #fig = plt.figure(figsize=(12,9))
    fig = plt.figure(figsize=(fig_width,fig_height))
    
    ax = plt.gca()
    #plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
    plt.title(r'171$\AA$  Visual: %s w/ Points' % (titles_vis[i]), y = 1.01, fontsize=25)  # no date / wavelength
    #im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
    im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
    plt.scatter(m2, l2, s=50, c='white')
    plt.xlim(0, vis[0].shape[1])
    plt.ylim(0, vis[0].shape[0])
    #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
    #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('Intensity', size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=17, pad=5) 
    #plt.tight_layout()
    #plt.savefig('%s/%s_%i_visual_%s.jpeg' % (path_name, date, wavelength, names_vis[i]))
    #plt.savefig('C:/Users/Brendan/Desktop/use_these/171_visual_%s_points.pdf' % names_vis[i], format='pdf')