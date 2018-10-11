# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 10:57:12 2017

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
#date = '20141227'

date = '20130626'    

h5 = np.load('%s/DATA/Output/%s/1700/param.npy' % (directory, date))
vis = np.load('%s/DATA/Output/%s/1700/visual.npy' % (directory, date))
vis = vis[1:-1,1:-1]


# trim x/y dimensions equally so that resulting region is 1600x1600    
trim_y = (h5.shape[1]-1600)//2
trim_x = (h5.shape[2]-1600)//2
h5 = h5[:, trim_y:h5.shape[1]-trim_y, trim_x:h5.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)  

trim_yv = (vis.shape[0]-1600)//2
trim_xv = (vis.shape[1]-1600)//2
vis = vis[trim_yv:vis.shape[0]-trim_yv, trim_xv:vis.shape[1]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

#date = '20140910'
#path_name = 'F:/Users/Brendan/Desktop/SolarProject/data/20130626'
path_name = 'C:/Users/Brendan/Desktop/same_scale'

# create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
titles = [r'Power Law Slope-Coefficient -- [$A$]', r'Power Law Index -- [$n$]', r'Power Law Tail -- [$C$]', r'Gaussian Amplitude -- [$\alpha$]', r'Gaussian Location [Seconds] -- [$\beta$]', r'Gaussian Width -- [$\sigma$]', 'F-Statistic', r'Gaussian Amplitude Scaled -- [$\alpha$]', 'P-Value']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']

heatmap = [h5]
wavelengths = [1700]

year = date[0:4]
month = date[4:6]
day = date[6:8]
date_title = '%s-%s-%s' % (year,month,day)


vmin = np.zeros((7))
vmax = np.zeros((7))


fig = plt.figure(figsize=(12,10))
ax1 = plt.gca()
ax1 = plt.subplot2grid((11,11),(0, 6), colspan=5, rowspan=5)


directory = 'D:'
date = '20130626'
wavelength = 1700
    
# create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'Power Law Index - $n$', r'Power Law Tail - $C$', r'Gaussian Amplitude [flux] - $\alpha$', r'Gaussian Location [min] - $\beta$', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', 'p-Value']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [min]', 'Width', 'F-Statistic', 'Amplitude Scaled', 'p-Value']

# load parameter array and visual images from file tree structure 
heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
visual = np.load('%s/DATA/Output/%s/%i/visual.npy'% (directory, date, wavelength))  
h_map = heatmaps

font_size = 19  # set the font size to be used for all text - titles, tick marks, text, labels

wavelength = wavelength    

visual = visual[1:-1,1:-1]
trim_yv = (visual.shape[0]-1600)//2
trim_xv = (visual.shape[1]-1600)//2
visual = visual[trim_yv:visual.shape[0]-trim_yv, trim_xv:visual.shape[1]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)

#"""
# trim x/y dimensions equally so that resulting region is 1600x1600    
trim_y = (h_map.shape[1]-1600)//2
trim_x = (h_map.shape[2]-1600)//2
h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

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


plt.rcParams["font.family"] = "Times New Roman"
font_size = 12  # set the font size to be used for all text - titles, tick marks, text, labels

# generate visual images
titles_vis = ['Average', 'Middle-File']
names_vis = ['average', 'mid']

vis = visual
#vis = vis[:,145:475,510:925]
#v_mask = v_mask[145:475,510:925]
#p_mask = p_mask[145:475,510:925]
#p_maskI = p_maskI[145:475,510:925]

#m1 = 200
#m2 = 410
#n1 = 600
#n2 = 822
m1 = 200
m2 = 410
n1 = 680
n2 = 902

#"""
#x_ticks = [0,100,200]
x_ticks = [20,120,220]
y_ticks = [0,100,200]  
#x_ind = [-200,-100,0]
x_ind = [-100,0,100]
y_ind = [-600,-500,-400]    
#"""

vis = vis[m1:m2,n1:n2]
v_mask = v_mask[m1:m2,n1:n2]
p_mask = p_mask[m1:m2,n1:n2]
p_maskI = p_maskI[m1:m2,n1:n2]
loc_mask = loc_mask[m1:m2,n1:n2]
#trim_yv = (vis.shape[1]-1600)/2
#trim_xv = (vis.shape[2]-1600)/2
#vis = vis[:, trim_yv:vis.shape[1]-trim_yv, trim_xv:vis.shape[2]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

#vis = vis[:, 0:vis.shape[1]-50, 0:500]  # for 20130626 blobs     

i = 0
v_min = np.percentile(vis,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
v_max = np.percentile(vis,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)     

delta = 1.
x = np.arange(0., v_mask.shape[1], delta)
y = np.arange(0., v_mask.shape[0], delta)
X, Y = np.meshgrid(x, y)
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
            
#ax1 = plt.gca()
#ax1 = plt.subplot2grid((2,11),(0, 6), colspan=5, rowspan=1)
#plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
plt.title('b) Visual %s' % (titles_vis[i]), y = 1.01, fontsize=font_size)  # no date / wavelength
#im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
#im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#cmap = cm.get_cmap('sdoaia%i' % wavelength, 10)    
#im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#im = ax.imshow(vis[i], cmap='gray', vmin = v_min, vmax = v_max)
im = ax1.imshow(vis, cmap='gray', vmin = 0.2, vmax = 2.2)
#CS = plt.contour(X, Y, Z, levels=[2.], linewidths=2, colors='white', linestyles='solid')
#CS = plt.contour(X, Y, Z, levels=[2.], linewidths=2, colors='black', linestyles='dashed')
CS = plt.contour(X, Y, Z, levels=[2.], linewidths=2, colors='red', linestyles='solid')
#im = ax.imshow(np.flipud(v_mask), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)

plt.xticks(x_ticks,x_ind,fontsize=font_size)
plt.yticks(y_ticks,y_ind,fontsize=font_size)
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


#h_min = np.percentile((1./np.exp(h_map[4]))/60.,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
#h_max = np.percentile((1./np.exp(h_map[4]))/60.,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
h_min = 3.
h_max = 5.
cmap = cm.get_cmap('jet_r', 10)

#h_range = np.abs(h_max-h_min)
h_range = np.abs(5.-3.)
h_step = h_range / 10.
c_ticks = np.zeros((11))
for h in range(11):
    c_ticks[h] = h_min + h_step*h 

#fig = plt.figure(figsize=(12,10))
#fig = plt.figure(figsize=(fig_width,fig_height))

ax2 = plt.gca()
ax2 = plt.subplot2grid((11,11),(0, 0), colspan=5, rowspan=5)
#plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
plt.title('a) Lorentzian Location', y = 1.01, fontsize=font_size)  # no date / wavelength
#im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
#im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#cmap = cm.get_cmap('sdoaia%i' % wavelength, 10)    
#im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#im = ax.imshow(loc_mask, cmap=cmap, vmin = h_min, vmax = h_max)
im = ax2.imshow(loc_mask, cmap=cmap, vmin = h_min, vmax = h_max)
#im = ax.imshow(loc_mask, cmap='gray', vmin = h_min, vmax = h_max)
#im = ax.imshow(np.flipud(v_mask), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
plt.xticks(x_ticks,x_ind,fontsize=font_size)
plt.yticks(y_ticks,y_ind,fontsize=font_size)
#ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax, format='%0.1f')
cbar.set_label(label='Lorentz. Loc. [min]', fontsize=font_size)
#cbar.set_label('Intensity', size=20, labelpad=10)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
cbar.set_ticks(c_ticks)

#CS = ax.contour(X1, Y1, Z1, levels=[2.], linewidths=3, colors='red', linestyles='solid')

#plt.tight_layout()
#plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
#plt.savefig('C:/Users/Brendan/Desktop/20130626_1700_gauss_loc.pdf', format='pdf', bbox_inches='tight')


    

plt.rc('mathtext', default='regular')
plt.rc('image', origin='lower', interpolation='nearest', cmap='gray')


directory = 'D:'
date = '20130626'
f_name = 'hmi_lev1_continuum_2013_06_26t06_00_00_45z_image_lev1'

image = '%s/FITS/%s/continuum/%s.fits' % (directory, date, f_name)


m1 = Map('%s' % image)



#v_mask = v_mask[5:-10,1:]
v_mask = v_mask[:-10,:]

###########
# plate scaling

#continuum_scale_factor = m1.scale.x / (0.609373 * u.arcsec)  #1600
continuum_scale_factor = m1.scale[0].value / (0.612898 * u.arcsec)  #1700

hmimap = m1.rotate(scale=continuum_scale_factor.value)

h1 = hmimap.data

#h1 = h1[1448:1658,1950:2162]
h1 = h1[1448:1658,1940:2162]  #revised 10/20


#############
############
###########

h_range = np.abs(60. - 10.)
h_step = h_range / 10.
c_ticks = np.zeros((11))
for h in range(11):
    c_ticks[h] = 10. + h_step*h 

print(np.percentile(h1, 1), np.percentile(h1, 99))   
# Load an example image
img = h1
img /= img.max()

# Contrast stretching
p2 = np.percentile(img, 1)
p98 = np.percentile(img, 99)
img_rescale = exposure.rescale_intensity(img, in_range=(p2, p98))

# Equalization
img_eq = exposure.equalize_hist(img)

# Adaptive Equalization
img_adapteq = exposure.equalize_adapthist(img, clip_limit=0.03)

delta = 1.
x = np.arange(0., img.shape[1], delta)
y = np.arange(0., img.shape[0], delta)
X, Y = np.meshgrid(x, y)

Z = img

mask = np.copy(img_adapteq)    
mask[img > 0.75] = np.NaN

 
#plt.figure(figsize=(12,10))
ax3 = plt.gca()
ax3 = plt.subplot2grid((11,11),(6, 6), colspan=5, rowspan=5)
#ax = plt.gca()  # get current axis -- to set colorbar 
#plt.title('2013/06/26: HMI Continuum | p-Value Mask Overlay', y=1.01, fontsize=19)
plt.title('d) HMI Continuum', y=1.01, fontsize=font_size)

im = ax3.imshow((h1/1e3),origin='lower', vmin=10., vmax=60)
img0 = img_as_float(img_rescale)
ax3.imshow(img0, cmap=plt.cm.gray)
img = img_as_float(mask)
ax3.imshow(mask, cmap=plt.cm.gray)

plt.xticks(x_ticks,x_ind,fontsize=font_size)
plt.yticks(y_ticks,y_ind,fontsize=font_size)
divider = make_axes_locatable(ax3)  # set colorbar to heatmap axis
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax, format='%i')
cbar.set_label(label='$I_{\mathrm{c}}$ [kDN/s]', fontsize=font_size)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
cbar.set_ticks(c_ticks)

x1 = np.arange(0., v_mask.shape[1], delta)
y1 = np.arange(0., v_mask.shape[0], delta)
X1, Y1 = np.meshgrid(x1, y1)
Z1 = v_mask

CS = ax3.contour(X1, Y1, Z1, levels=[2.], linewidths=2, colors='red', linestyles='solid')

#plt.savefig('C:/Users/Brendan/Desktop/20130626_1700_continuum_pvalue_overlay_scaled_to_aia_equalization.pdf', format='pdf', bbox_inches='tight')




directory = 'D:'
date = '20130626'
f_name = 'hmi_lev1_magnetogram_2013_06_26t06_00_00_45z_image_lev1'

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
h1 = h1[1448:1658,1940:2162]


hmimag = plt.get_cmap('hmimag')

h_range = np.abs(2.)
h_step = h_range / 10.
c_ticks = np.zeros((11))
for h in range(11):
    c_ticks[h] = -1. + h_step*h 


#fig = plt.figure(figsize=(12,10))
#ax = fig.add_subplot(111)
ax4 = plt.gca()
ax4 = plt.subplot2grid((11,11),(6, 0), colspan=5, rowspan=5)
#plt.title('2013/06/26: HMI Magnetogram | p-Value Mask Overlay', y=1.01, fontsize=19)
plt.title('c) HMI Magnetogram', y=1.01, fontsize=font_size)
#im = plt.imshow((h1/1e3), cmap=hmimag,vmin=-3,vmax=3)
#im = plt.imshow(np.fliplr(np.flipud(d2)),origin='lower',cmap=hmimag,vmin=-3000,vmax=3000)
#im = ax4.imshow((h1/1e3), cmap='gray',vmin=-1,vmax=1)
im = ax4.imshow((h1/1e3), cmap='gray',vmin=-1.,vmax=1.)
plt.xticks(x_ticks,x_ind,fontsize=font_size)
plt.yticks(y_ticks,y_ind,fontsize=font_size)

divider = make_axes_locatable(ax4)  # set colorbar to heatmap axis
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax, format='%0.1f')
cbar.set_label(label='$B_{\mathrm{los}}$ [kG]', fontsize=font_size)
cbar.ax.tick_params(labelsize=font_size, pad=5) 
cbar.set_ticks(c_ticks)

x1 = np.arange(0., v_mask.shape[1], delta)
y1 = np.arange(0., v_mask.shape[0], delta)
X1, Y1 = np.meshgrid(x1, y1)
#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
#Z = 10.0 * (Z2 - Z1)
Z1 = v_mask

CS = ax4.contour(X1, Y1, Z1, levels=[2.,], linewidths=2, colors='red', linestyles='solid')
#CS = ax.contour(X1, Y1, Z1, levels=[325.], linewidths=5, colors='blue', linestyles='solid')
#CS = ax.contourf(X1, Y1, Z1, levels=[2.,100000.], colors='red', hatches=['/'], extend='lower')
#CS = ax.contour(X1, Y1, Z1, levels=[2.], linewidths=3, colors='black', linestyles='dashed')
plt.savefig('C:/Users/Brendan/Desktop/20130626_1700_2x2_figureF.pdf', format='pdf', bbox_inches='tight')
#"""