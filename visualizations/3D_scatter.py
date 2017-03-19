# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 18:33:16 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

#HEATMAPS = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20130530_1600_2300_2600_2200_3000_float_numpy.npy')
HEATMAPS = np.load('C:/Users/Brendan/Desktop/SDO/param_20130530_1600_2300_2600i_2200_3000j_data_rebin4b.npy')

titles = ['Power Law Slope Coefficient', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width', '($/chi^2$)']
names = ['PL_A', 'Slopes', 'PL_C', 'Gauss_Amp', 'Gauss_Loc', 'Gauss_Wid', 'Chi2']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '($/chi^2$)']
#vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore
#vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore
#wavelength = wavelength
#year = date[0:4]
#month = date[4:6]
#day = date[6:8]
wavelength = 1600
year = '2013'
month = '05'
day = '30'
date_title = '%s-%s-%s' % (year,month,day)

h_map = HEATMAPS
h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]

flat_hmap = np.zeros((h_map.shape[0],h_map.shape[1]*h_map.shape[2]))


for i in range(h_map.shape[0]):
    if i == 6:
        NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
        flat_hmap[i] = np.reshape(NaN_replace, (h_map.shape[1]*h_map.shape[2]))
    else:
        flat_hmap[i] = np.reshape(h_map[i], (h_map.shape[1]*h_map.shape[2]))


fig = plt.figure(figsize=(15,9))
ax = fig.add_subplot(111,projection='3d')
#ax = fig.add_subplot(111)
#ax.scatter(fl_plC, fl_slopes, fl_plA,marker='.')
i = 0
j = 1
k = 2
xmin = np.percentile(flat_hmap[i],1)
xmax = np.percentile(flat_hmap[i],99)
ymin = np.percentile(flat_hmap[j],1)
ymax = np.percentile(flat_hmap[j],99)
zmin = np.percentile(flat_hmap[k],1)
zmax = np.percentile(flat_hmap[k],99)
x_margin = (xmax-xmin)*0.05
y_margin = (ymax-ymin)*0.05
z_margin = (zmax-zmin)*0.05
ax.set_xlim(xmin-x_margin, xmax+x_margin)
ax.set_ylim(ymin-y_margin, ymax+y_margin)
ax.set_zlim(zmin-z_margin, zmax+z_margin)
ax.set_xlabel('%s' % titles[i])
ax.set_ylabel('%s' % titles[j])
ax.set_zlabel('%s' % titles[k])
ax.scatter(flat_hmap[i], flat_hmap[j], flat_hmap[k], marker='.')            
            
"""
for j in range(h_map.shape[0]):
    for i in range(h_map.shape[0]):    
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
            path_name='C:/Users/Brendan/Desktop/PHYS 326/test_temp'
            date = '20130530'
            plt.savefig('%s/%s_%i_scatter_%s_vs_%s.jpeg' % (path_name, date, wavelength, names[i], names[j]))
"""