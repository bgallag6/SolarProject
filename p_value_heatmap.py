# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 10:41:20 2017

@author: Brendan
"""
from timeit import default_timer as timer
from scipy.stats import f
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import distributions
from mpl_toolkits.axes_grid1 import make_axes_locatable

#h_map = np.load('C:/Users/Brendan/Desktop/SDO/20130815_193_1000_1600i_1950_2950j_rebin2_params_mpi.npy')
h_map = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_171_-500_500i_-500_600j_param_slope6_arthm.npy')
  
df1, df2 = 3, 6
x = distributions.f.ppf(0.95,df1,df2)
flat_param = np.reshape(h_map[6], (h_map[6].shape[0]*h_map[6].shape[1]))

p_value = f.cdf(h_map[6], df1, df2)


y = f.sf(h_map[6], df1, df2) # equals 1-p_value

yy = f.sf(15, df1, df2)


#mask_arr = np.copy(y)

p_arr = [0.15, 0.1, 0.05, 0.01, 0.005, 0.001]

"""
fig = plt.figure(figsize=(13,9))
ax1 = plt.subplot2grid((1,11),(0, 0), colspan=5, rowspan=1)
ax1 = plt.gca()  # get current axis -- to set colorbar 
ax1.set_title('2013/08/15 193A, F-Value', y = 1.01, fontsize=17)
#plt.title('2013/08/15 193A, F-Value', y = 1.01, fontsize=20)
im = ax1.imshow(h_map[6])
divider = make_axes_locatable(ax1)  # set colorbar to heatmap axis
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
#cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
cbar.ax.tick_params(labelsize=17, pad=5) 
"""
start = timer()

for n in range(len(p_arr)):
    c = p_arr[n]

    mask_arr = np.copy(y)
    
    mask_arr[y > c] = np.NaN  # this replaces the loop below!  about 5x faster

    #for i in range(y.shape[0]):
    #    for j in range(y.shape[1]):
    #        if y[i][j] > c:
    #            mask_arr[i][j] = np.NaN
    #"""            
    fig = plt.figure(figsize=(16,6))
    ax1 = plt.subplot2grid((1,11),(0, 0), colspan=5, rowspan=1)
    ax1 = plt.gca()  # get current axis -- to set colorbar 
    ax1.set_title('2013/08/15 193A, F-Value', y = 1.01, fontsize=15)
    #plt.title('2013/08/15 193A, F-Value', y = 1.01, fontsize=20)
    im = ax1.imshow(h_map[6])
    divider = make_axes_locatable(ax1)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=15, pad=5) 

    #fig = plt.figure(figsize=(13,9))
    ax2 = plt.subplot2grid((1,11),(0, 6), colspan=5, rowspan=1)
    ax2 = plt.gca()  # get current axis -- to set colorbar 
    #plt.title('2013/08/15 193A, P-Value [Mask > %f]' % c, y = 1.01, fontsize=20)
    ax2.set_title('2013/08/15 193A, P-Value [Mask > %f]' % c, y = 1.01, fontsize=15)
    im = ax2.imshow(mask_arr)
    divider = make_axes_locatable(ax2)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax, format='%.3f')
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=15, pad=5) 
    plt.tight_layout()
    #plt.savefig('C:/Users/Brendan/Desktop/20130815_193_pvalue_%.3f.jpeg' % c)
    #"""

   
cf = 0.005
mask_arr = np.copy(y)
mask_amp = np.copy(h_map[3])
mask_loc = np.copy(h_map[4])
mask_wid = np.copy(h_map[5])
"""
for i in range(y.shape[0]):
        for j in range(y.shape[1]):
            if y[i][j] > cf:
                mask_arr[i][j] = np.NaN
                mask_amp[i][j] = np.NaN
                mask_loc[i][j] = np.NaN
                mask_wid[i][j] = np.NaN
"""
mask_arr[mask_arr > cf] = np.NaN
mask_amp[mask_arr > cf] = np.NaN
mask_loc[mask_arr > cf] = np.NaN
mask_wid[mask_arr > cf] = np.NaN
                
"""
fig = plt.figure(figsize=(10,7))
plt.imshow(mask_arr)
fig = plt.figure(figsize=(10,7))
plt.imshow(h_map[3])
fig = plt.figure(figsize=(10,7))
plt.imshow(mask_amp)
fig = plt.figure(figsize=(10,7))
plt.imshow(h_map[4])
fig = plt.figure(figsize=(10,7))
plt.imshow(mask_loc)
fig = plt.figure(figsize=(10,7))
plt.imshow(h_map[5])            
fig = plt.figure(figsize=(10,7))
plt.imshow(mask_wid)     
"""          

T_final = timer() - start
T_min_final, T_sec_final = divmod(T_final, 60)
T_hr_final, T_min_final = divmod(T_min_final, 60)
print "Total program time = %f sec" % T_final 
    

