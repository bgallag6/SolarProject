# -*- coding: utf-8 -*-
"""
Created on Sun Mar 25 21:31:52 2018

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm
from scipy import stats
import matplotlib
import scipy

def Gauss(f, P, fp, fw):
    #return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2) + C
    return P*np.exp(-0.5*((f-fp)/fw)**2)


directory = 'F:'
date = '20130626'
wavelength = 171

#directory = 'S:'
#date = '20001111'
#wavelength = 1600

matplotlib.rc('text', usetex = True)  # use with latex commands

# 11-param-list
titles = [r'Power Law Slope-Coefficient [flux] - A', r'(b) Power Law Index n', r'Power Law Tail - C', r'Gaussian Amplitude [flux] - α', r'(c) Gauss. Loc. β [min]', r'Gaussian Width - σ', 'F-Statistic', r'Gaussian Amplitude Scaled - α', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
#titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'(b) Power Law Index $n$', r'Power Law Tail - C', r'Gaussian Amplitude [flux] - $\alpha$', r'(c) Gauss. Loc. $\beta$ [min]', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [min]', 'Width', 'F-Statistic', 'Amplitude Scaled', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
names = ['slope_coeff', 'index', 'tail', 'lorentz_amp', 'lorentz_loc', 'lorentz_wid', 'f_test', 'lorentz_amp_scaled', 'r_value', 'roll_freq', 'chisqr']

# load parameter array and visual images from file tree structure 
heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
visual = np.load('%s/DATA/Output/%s/%i/visual.npy'% (directory, date, wavelength))  

#visual = visual[1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
visual = visual[1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
h_map = heatmaps    

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels
plt.rc('font',**{'family':'serif','serif':['Times']})

wavelength = wavelength    

#"""
# trim x/y dimensions equally so that resulting region is 1600x1600    
trim_y = int((h_map.shape[1]-1600)/2)
trim_x = int((h_map.shape[2]-1600)/2)
h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
x_ind = [-800,-600,-400,-200,0,200,400,600,800]
y_ind = [800,600,400,200,0,-200,-400,-600,-800]    
#"""

# generate p-value heatmap + masked Gaussian component heatmaps
df1, df2 = 3, 6  # degrees of freedom for model M1, M2
p_val = ff.sf(h_map[6], df1, df2)

#p_mask = np.copy(p_val)

mask_thresh = 0.005  # significance threshold - masked above this value
   
p_mask = np.copy(p_val)
amp_mask = np.copy(h_map[3])
loc_mask = np.copy(h_map[4])
wid_mask = np.copy(h_map[5])    

# mask the Gaussian component arrays with NaNs if above threshold 
p_mask[p_val > mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
amp_mask[p_val > mask_thresh] = np.NaN
loc_mask[p_val > mask_thresh] = np.NaN
wid_mask[p_val > mask_thresh] = np.NaN

"""
# for creating sunspot umbra + PPV contour overlays from 1600
v1600 = np.load('%s/DATA/Output/%s/1600/visual.npy'% (directory, date))  
v1600 = v1600[:,1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
p1600 = np.load('%s/DATA/Output/%s/1600/param.npy'% (directory, date)) 

h_map = h_map[:,:p1600.shape[1],:p1600.shape[2]]
visual = visual[:,:v1600.shape[1],:v1600.shape[2]]

#v_mask = np.copy(visual[0])
p1600_val = ff.sf(p1600[6], df1, df2)
#p1600_mask = np.copy(p1600_val)
v_mask = np.copy(v1600[0])
   
v_mask[p1600_val < mask_thresh] = 1.  # invert mask, set equal to 1. -- so can make contour
""" 


# determine percentage of region masked 
count = np.count_nonzero(np.isnan(p_mask))   
total_pix = p_val.shape[0]*p_val.shape[1]
mask_percent = ((np.float(count))/total_pix)*100
            
loc_mask = (1./np.exp(loc_mask))/60.  # convert Gaussian location to minutes
plots = [p_mask, amp_mask, loc_mask, wid_mask]  # make array of masked plots to iterate over


fig_width = 12
fig_height = 10


fig = plt.figure(figsize=(fig_width,fig_height))
ax = plt.gca()  # get current axis -- to set colorbar 
plt.title(r'(b) Power Law Index $n$', y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength

h_min = np.percentile(h_map[1],1)
h_max = np.percentile(h_map[1],99)
cmap = cm.get_cmap('jet', 10)

# specify colorbar ticks to be at boundaries of segments
h_range = np.abs(h_max-h_min)
h_step = h_range / 10.
c_ticks = np.zeros((11))
for h in range(11):
    c_ticks[h] = h_min + h_step*h 
    
im = ax.imshow(np.flipud(h_map[1]), cmap = cmap, vmin=h_min, vmax=h_max)
plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
cax = divider.append_axes("right", size="3%", pad=0.07)

cbar = plt.colorbar(im,cax=cax, format='%0.2f')
cbar.ax.tick_params(labelsize=font_size, pad=5) 
cbar.set_ticks(c_ticks)

#plt.savefig('C:/Users/Brendan/Desktop/%s_%i_index.pdf' % (date, wavelength), format='pdf', bbox_inches='tight')



#"""
flat_param = np.reshape(h_map[1], (h_map[1].shape[0]*h_map[1].shape[1]))

# calculate some statistics
mean = np.mean(flat_param)
sigma = np.std(flat_param)   
print('Index: mean = %0.2f, sigma = %0.2f' % (mean, sigma))

fig = plt.figure(figsize=(fig_width+1,fig_height))
plt.title('%i \AA\\ | Power Law Index' % wavelength, y = 1.02, fontsize=font_size)  # no date / wavelength
plt.xlabel('%s' % cbar_labels[1], fontsize=font_size, labelpad=10)
plt.ylabel('Bin Count', fontsize=font_size, labelpad=10)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
#plt.xlim(h_min, h_max)
y, x, _ = plt.hist(flat_param, bins=200, edgecolor='black')
#y, x, _ = plt.hist(flat_param, bins=100, range=(h_min, h_max), edgecolor='black')


n=y[1:-2]
bins=x[1:-2]
elem = np.argmax(n)
bin_max = bins[elem]
plt.ylim(0, y.max()*1.1)

f = x[:-1]
s = y
#f = x[:-26]  # for 1700 so will compute gaussian fit
#s = y[:-25]

p_low = [0.,-1.,0.]
p_high = [1e6, 2., 5.]
#nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(Gauss, f, s, method='dogbox', max_nfev=10000)     
nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(Gauss, f, s, bounds=(p_low, p_high))       
#P, fp, fw, C = nlfit_gp  # unpack fitting parameters
P, fp, fw = nlfit_gp  # unpack fitting parameters          
#g_fit = Gauss(f, P,fp,fw, C)  
g_fit = Gauss(f, P,fp,fw)       
gauss_center = np.exp(fp)
gauss_wid = fw

plt.plot(f,s, linewidth=1.5)
plt.xlim(h_min,h_max)
plt.plot(f,g_fit, linestyle='dashed', linewidth=2.)
plt.vlines(gauss_center,0,y.max()*1.1, linestyle='dashed', color='red', linewidth=2., label='center=%0.4f' % gauss_center)
plt.vlines(0, 0, y.max()*1.1, color='white', linestyle='dashed', linewidth=1.5, label='width=%0.4f' % gauss_wid)


plt.vlines(bin_max, 0, y.max()*1.1, color='black', linestyle='dotted', linewidth=2., label='mode=%0.4f' % bin_max)  
plt.vlines(0, 0, y.max()*1.1, color='white', linestyle='dashed', linewidth=1.5, label='sigma=%0.4f' % sigma)

legend = plt.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
for label in legend.get_lines():
    label.set_linewidth(2.0)

#plt.savefig('C:/Users/Brendan/Desktop/20130626_%i_index_histB.pdf' % wavelength, format='pdf', bbox_inches='tight')
   
    
## r values
r_mean = np.mean(h_map[8])
r_sigma = np.std(h_map[8])   
print('R-value: mean = %0.2f, sigma = %0.2f' % (r_mean, r_sigma))


## percentage of pixels using M1 over M2
m1count = np.count_nonzero(np.isnan(h_map[4]))   
total_pix = h_map[4].shape[0]*h_map[4].shape[1]
m1percent = ((np.float(m1count))/total_pix)*100
print('M1 percent = %0.2f' % m1percent)


print(np.min(h_map[0]))
print(np.min(h_map[1]))
print(np.min(h_map[2]))
