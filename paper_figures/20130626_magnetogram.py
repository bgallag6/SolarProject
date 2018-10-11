# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 18:50:01 2017

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
from matplotlib.colors import LogNorm
from skimage import data, img_as_float
from skimage import exposure


directory = 'D:'
date = '20130626'
wavelength = 1700
    
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

font_size = 19  # set the font size to be used for all text - titles, tick marks, text, labels

wavelength = wavelength    

visual = visual[1:-1,1:-1]
#h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]  # trim last row and column from array (originally needed since went one past)
trim_yv = int((visual.shape[0]-1600)/2)
trim_xv = int((visual.shape[1]-1600)/2)
visual = visual[trim_yv:visual.shape[0]-trim_yv-1, trim_xv:visual.shape[1]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)

#"""
# trim x/y dimensions equally so that resulting region is 1600x1600    
trim_y = int((h_map.shape[1]-1600)/2)
trim_x = int((h_map.shape[2]-1600)/2)
h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

"""
x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
x_ind = [-800,-600,-400,-200,0,200,400,600,800]
y_ind = [800,600,400,200,0,-200,-400,-600,-800]    
"""

#"""
x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
x_ind = [-800,-600,-400,-200,0,200,400,600,800]
y_ind = [800,600,400,200,0,-200,-400,-600,-800]  
y_ind = np.flipud(y_ind)  

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

v_mask = np.copy(visual)
#v_mask = v_mask[1:-1,1:-1]

# mask the Gaussian component arrays with NaNs if above threshold 
p_mask[p_val > mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
p_maskI[p_val < mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
v_mask[p_val < mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
v_mask[p_val < mask_thresh] = 1.  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN

loc_mask = np.copy(h_map[4])
loc_mask[p_val > mask_thresh] = np.NaN             
loc_mask = (1./np.exp(loc_mask))/60.  # convert Gaussian location to minutes

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


plt.rcParams["font.family"] = "Times New Roman"
font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels

# generate visual images
titles_vis = ['Average', 'Middle-File']
names_vis = ['average', 'mid']

vis = visual
#vis = vis[:,145:475,510:925]
#v_mask = v_mask[145:475,510:925]
#p_mask = p_mask[145:475,510:925]
#p_maskI = p_maskI[145:475,510:925]

m1 = 0
m2 = 1600
n1 = 0
n2 = 1600

#"""
x_ticks = [0,100,200]
y_ticks = [0,100,200]  
x_ind = [-200,-100,0]
y_ind = [-600,-500,-400]    
#"""

#vis = vis[:,m1:m2,n1:n2]

#trim_yv = (vis.shape[1]-1600)/2
#trim_xv = (vis.shape[2]-1600)/2
#vis = vis[:, trim_yv:vis.shape[1]-trim_yv, trim_xv:vis.shape[2]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

#vis = vis[:, 0:vis.shape[1]-50, 0:500]  # for 20130626 blobs     
x_ticks = [0,200,400,600,700,750,800,850,900,1000,1200,1400,1600]
y_ticks = [0,200,300,400,600,800,1000,1200,1400,1600]  
x_ind = [-800,-600,-400,-200,-100,-50,0,50,100,200,400,600,800]
y_ind = [800,600,400,200,0,-200,-400,-500,-600,-800]
y_ind = np.flipud(y_ind)

i = 0
v_min = np.percentile(vis,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
v_max = np.percentile(vis,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)     

delta = 1.
x = np.arange(0., v_mask.shape[1], delta)
y = np.arange(0., v_mask.shape[0], delta)
X, Y = np.meshgrid(x, y)
#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
#Z = 10.0 * (Z2 - Z1)
Z = v_mask



vis = vis/1e3
v_min = v_min/1e3
v_max = v_max/1e3

print(v_min, v_max)

h_range = np.abs(2.2 - 0.2)
h_step = h_range / 10.
c_ticks = np.zeros((11))
for h in range(11):
    c_ticks[h] = 0.2 + h_step*h 
          
          
fig = plt.figure(figsize=(12,10))  
ax1 = plt.gca()       
plt.title('b) Visual %s' % (titles_vis[i]), y = 1.01, fontsize=font_size)  # no date / wavelength
im = ax1.imshow(vis, cmap='sdoaia1700', vmin = 0.2, vmax = 2.2)
plt.xticks(x_ticks,x_ind,fontsize=font_size)
plt.yticks(y_ticks,y_ind,fontsize=font_size)
#plt.ylabel('Y-position [pixels]', fontsize=15)
#plt.xlabel('X-position [pixels]', fontsize=15)
#ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="3%", pad=0.07)
#cbar = plt.colorbar(im,cax=cax)
cbar = plt.colorbar(im,cax=cax, format='%0.1f')
cbar.set_label(label='Intensity [kDN/s]', fontsize=font_size)
#cbar.set_label('Intensity', size=20, labelpad=10)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
cbar.set_ticks(c_ticks)
#plt.tight_layout()
#plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
#plt.savefig('C:/Users/Brendan/Desktop/20130626_1700_visual_average_pvalue_mask_contour_overlay.pdf', format='pdf', bbox_inches='tight')
#"""


plt.rc('mathtext', default='regular')
plt.rc('image', origin='lower', interpolation='nearest', cmap='gray')

directory = 'D:'
date = '20130626'
f_name = 'hmi_lev1_magnetogram_2013_06_26t06_00_00_45z_image_lev1'
#f_name = 'hmi_lev1_magnetogram_2013_06_26t00_00_00_45z_image_lev1'

image = '%s/FITS/%s/magnetogram/%s.fits' % (directory, date, f_name)

m1 = Map('%s' % image)

###########
# plate scaling


#m1 = 145 -- 200
#m2 = 475 -- 410
#n1 = 510 -- 610
#n2 = 925 -- 822

#hmi_scale_factor = m1.scale.x / (0.6 * u.arcsec)
#hmi_scale_factor = m1.scale.x / (0.609373 * u.arcsec)  #1600
hmi_scale_factor = m1.scale[0].value / (0.612898 * u.arcsec)  #1700

hmimap = m1.rotate(scale=hmi_scale_factor.value)

h1 = hmimap.data

#h1 = h1[1452:1662,1951:2163]
#h1 = h1[1250:2850,1250:2850]  # for initial image
h1 = h1[1248:2848,1260:2860]  # for middle image



hmimag = plt.get_cmap('hmimag')

h_range = np.abs(2.)
h_step = h_range / 10.
c_ticks = np.zeros((11))
for h in range(11):
    c_ticks[h] = -1. + h_step*h 


fig = plt.figure(figsize=(12,10))  
ax4 = plt.gca()
plt.title('(b) HMI Magnetogram', y=1.02, fontsize=font_size, fontname="Times New Roman")
#im = ax4.imshow((h1/1e3), cmap='hmimag',vmin=-1.,vmax=1.)
im = ax4.imshow((h1/1e3), cmap='gray',vmin=-1.,vmax=1.)
plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
ax4.tick_params(axis='both', which='major', pad=10)
#plt.ylabel('Y-position [pixels]', fontsize=15)
#plt.xlabel('X-position [pixels]', fontsize=15)
divider = make_axes_locatable(ax4)  # set colorbar to heatmap axis
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax, format='%0.1f')
cbar.set_label(label='$B_{\mathrm{los}}$ [kG]', fontsize=font_size, labelpad=10)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
cbar.set_ticks(c_ticks)

#plt.savefig('C:/Users/Brendan/Desktop/20130626_magnetogram_fullF.pdf', format='pdf', bbox_inches='tight')
#"""