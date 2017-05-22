# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 18:08:06 2017

@author: Brendan
"""

import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy
from sunpy.map import Map
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm

directory = 'F:/Users/Brendan/Desktop/SolarProject'
#date = '20141227'

dates = ['20101208','20111210','20121018', '20131118', '20140112', '20140606', '20140818', '20140910', '20141227', '20150104','20160414', '20160426', '20160520', '20160905', '20170329']
contours = [85.,125., 95., 110., 67., 95., 75., 115., 75., 113., 87., 60., 60., 95., 83.]
circ1_center_x = [154,0, 0, 520.,0,  163.5, 99.,0,  154., 0, 540., 141., 327., 661., 512.]
circ2_center_x = [154,0, 0, 520.,0,  0, 99.,0,  154., 0, 0., 0., 327., 0., 0.]
circ1_center_y = [125,0, 0, 290.,0,  158., 81.,0,  131., 0, 433., 113., 245., 285., 350.]
circ2_center_y = [125,0, 0, 290.,0,  0, 81.,0,  131.,0,  0., 0., 245., 0., 0.]
circ1_radi = [40,0, 0, 55.,0,  50., 20.,0,  24., 0, 85., 26., 50., 51., 37.]
circ2_radi = [55,0, 0, 73.,0,  0., 39.,0,  38., 0, 0., 0., 73., 0., 0.]

l = 6
date = dates[l]    

#h1 = np.load('%s/DATA/Output/%s/171/param.npy' % (directory, date))
#h2 = np.load('%s/DATA/Output/%s/193/param.npy' % (directory, date))
#h3 = np.load('%s/DATA/Output/%s/211/param.npy' % (directory, date))
#h4 = np.load('%s/DATA/Output/%s/304/param.npy' % (directory, date))
h5 = np.load('%s/DATA/Output/%s/1600/param.npy' % (directory, date))
vis = np.load('%s/DATA/Output/%s/1600/visual.npy' % (directory, date))
vis = vis[:,1:-1,1:-1]

#date = '20140910'
#path_name = 'F:/Users/Brendan/Desktop/SolarProject/data/20130626'
path_name = 'C:/Users/Brendan/Desktop/same_scale'

# create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
titles = [r'Power Law Slope-Coefficient -- [$A$]', r'Power Law Index -- [$n$]', r'Power Law Tail -- [$C$]', r'Gaussian Amplitude -- [$\alpha$]', r'Gaussian Location [Seconds] -- [$\beta$]', r'Gaussian Width -- [$\sigma$]', 'F-Statistic', r'Gaussian Amplitude Scaled -- [$\alpha$]', 'P-Value']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']

#heatmap = [h1,h2,h3,h4,h5]
heatmap = [h5]
#wavelengths = [17100,1930,211,30400,160000]
wavelengths = [160000]
#wavelengths = [171,193,211,304,1600]

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
#for c in range(1):
    h_map = heatmap[c]
    wavelength = wavelengths[c]

    #h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]  # trim last row and column from array (originally needed since went one past)
    
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
    
   
    visual = vis[0]
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
    
    delta = 1.
    x = np.arange(0., visual.shape[1], delta)
    y = np.arange(0., visual.shape[0], delta)
    X, Y = np.meshgrid(x, y)
    #Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    #Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    # difference of Gaussians
    #Z = 10.0 * (Z2 - Z1)
    Z = visual
    
    v_min = np.percentile(visual,1)
    v_max = np.percentile(visual,99)
    
    #for i in range(7):
    for i in range(5,6):
        
        
        if i == 6:
            NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
            h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #cmap = 'jet'
            #cmap = cm.get_cmap('jet', 10)
            cmap = cm.get_cmap('jet')
        elif i == 4:
            h_map[i] = 1./(np.exp(h_map[i]))
            #h_map2[i] = 1./(np.exp(h_map2[i]))
            h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #cmap = 'jet_r'  # reverse color-scale for Gaussian Location, because of flipped frequencies to seconds
            #cmap = cm.get_cmap('jet_r', 10)
            cmap = cm.get_cmap('jet_r')
        elif i == 8:
            df1, df2 = 3, 6
            h_map[6] = ff.sf(h_map[6], df1, df2)
            h_min = np.percentile(h_map[6],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[6],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #cmap = 'jet'      
            #cmap = cm.get_cmap('jet', 10)      
            cmap = cm.get_cmap('jet')
        else:
            h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #cmap = 'jet'
            #cmap = cm.get_cmap('jet', 10)
            cmap = cm.get_cmap('jet')
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
        
        #seg = ['9_11', '11_13', '13_15', '15_17', '4hr_13_17', '8hr_9_17']        
        
        lw = 2.        
        
        #fig = plt.figure(figsize=(13,9))
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
        #plt.title('%s -- Segment %i' % (titles[i], c), y = 1.01, fontsize=25)  # no date / wavelength
        #plt.title('%s -- Segment Hours: %s' % (titles[i], seg[c]), y = 1.01, fontsize=25)  # no date / wavelength
        circle = plt.Circle((circ1_center_x[l], circ1_center_y[l]), circ1_radi[l], linewidth=lw, color='white', fill=False, linestyle='dashed')
        circleB = plt.Circle((circ1_center_x[l], circ1_center_y[l]), circ1_radi[l], linewidth=lw+2.5, color='black', fill=False, linestyle='dashed')
        circle2 = plt.Circle((circ2_center_x[l], circ2_center_y[l]), circ2_radi[l], linewidth=lw, color='white', fill=False, linestyle='dashed')
        circle2B = plt.Circle((circ2_center_x[l], circ2_center_y[l]), circ2_radi[l], linewidth=lw+2.5, color='black', fill=False, linestyle='dashed')
        ax.add_artist(circleB)
        ax.add_artist(circle)
        ax.add_artist(circle2B)
        ax.add_artist(circle2)
        CS = plt.contour(X, Y, Z, levels=[contours[l]], linewidths=lw, colors='black', linestyles='solid')
        CS = plt.contour(X, Y, Z, levels=[contours[l]], linewidths=lw, colors='white', linestyles='dashed')
        #im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=vmin[i], vmax=vmax[i])
        im = ax.imshow(h_map[i], cmap = cmap, vmin=vmin[i], vmax=vmax[i])
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
        #plt.savefig('%s/%s_%i_%s_same_%i.jpeg' % (path_name, date, wavelength, names[i], c))
        #plt.savefig('%s/%s_%s_same_%i.jpeg' % (path_name, date, names[i], wavelength))
        #plt.savefig('%s/%s_%s_same_%iA.pdf' % (path_name, date, names[i], wavelength), format='pdf')
        #plt.savefig('%s/%s_%i_%s_same_%s.pdf' % (path_name, date, wavelength, names[i], seg[c]), format='pdf')
        #plt.close()

        #plt.title(r'%s %i$\AA$: Visual Average' % (dates[i],wavelength) + '\n Umbra [solid] + Penumbra [dashed]', y=1.01, fontsize=21)

#fig = plt.figure(figsize=(13,9))
fig = plt.figure(figsize=(fig_width,fig_height))
v_min = np.percentile(visual,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
v_max = np.percentile(visual,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)   
ax = plt.gca()  # get current axis -- to set colorbar 
plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
#plt.title('%s -- Segment %i' % (titles[i], c), y = 1.01, fontsize=25)  # no date / wavelength
#plt.title('%s -- Segment Hours: %s' % (titles[i], seg[c]), y = 1.01, fontsize=25)  # no date / wavelength
circle = plt.Circle((circ1_center_x[l], circ1_center_y[l]), circ1_radi[l], linewidth=lw, color='white', fill=False, linestyle='dashed')
circleB = plt.Circle((circ1_center_x[l], circ1_center_y[l]), circ1_radi[l], linewidth=lw+2.5, color='black', fill=False, linestyle='dashed')
circle2 = plt.Circle((circ2_center_x[l], circ2_center_y[l]), circ2_radi[l], linewidth=lw, color='white', fill=False, linestyle='dashed')
circle2B = plt.Circle((circ2_center_x[l], circ2_center_y[l]), circ2_radi[l], linewidth=lw+2.5, color='black', fill=False, linestyle='dashed')
ax.add_artist(circleB)
ax.add_artist(circle)
ax.add_artist(circle2B)
ax.add_artist(circle2)
CS = plt.contour(X, Y, Z, levels=[contours[l]], linewidths=lw, colors='black', linestyles='solid')
CS = plt.contour(X, Y, Z, levels=[contours[l]], linewidths=lw, colors='white', linestyles='dashed')
#im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=vmin[i], vmax=vmax[i])
im = ax.imshow(visual, cmap = 'sdoaia1600', vmin=v_min, vmax=v_max)
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
#plt.savefig('%s/%s_%i_%s_same_%i.jpeg' % (path_name, date, wavelength, names[i], c))
#plt.savefig('%s/%s_visual_average_%i.jpeg' % (path_name, date, wavelength))
#plt.savefig('%s/%s_%s_same_%iB.pdf' % (path_name, date, names[i], wavelength), format='pdf')
#plt.savefig('%s/%s_%i_%s_same_%s.pdf' % (path_name, date, wavelength, names[i], seg[c]), format='pdf')
#plt.close()

#plt.title(r'%s %i$\AA$: Visual Average' % (dates[i],wavelength) + '\n Umbra [solid] + Penumbra [dashed]', y=1.01, fontsize=21)