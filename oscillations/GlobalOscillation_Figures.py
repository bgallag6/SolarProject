# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 16:39:55 2017

@author: Brendan
"""

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import scipy
from scipy import fftpack  # doesnt work in module when called here???
from matplotlib import cm
from numpy.random import randn
import matplotlib.colors as colors
from timeit import default_timer as timer
from scipy.stats import f
import matplotlib.patches as patches
from scipy.stats import f as ff
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sunpy
from sunpy.map import Map
import astropy.units as u

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    #return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)
    return A2*f2**-n2 + C2 + P2*(1./ ((np.pi*fw2)*(1.+((np.log(f2)-fp2)/fw2)**2)))
    

#directory = 'F:/Users/Brendan/Desktop/SolarProject'
directory = 'S:'
date = '20120111'
wavelength = 1700

#"""
DATA = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))

Ex =  np.load('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength))

TIME = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))
  
ex_interp = np.interp(TIME,Ex,Ex)  # interpolate pixel-intensity values onto specified time grid
#"""  

pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
print DATA.shape
       
#v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27        

#x0 = 800 # 20120111
#y0 = 800
x0 = 796 # 20120111
y0 = 797

#x0 = 970 # 20170111
#y0 = 970


timeseries = DATA[:,y0,x0] / ex_interp


### Plot selected pixel's time series

fig = plt.figure(figsize=(12,10))
ax = plt.gca()
plt.title(r'(b) Time Series', y = 1.01, fontsize=25)
ax.tick_params(axis='both', which='major', pad=10)
plt.xlim(0,240)
plt.plot(TIME/60., timeseries, 'k')
plt.xlabel(r'Time [min]', fontsize=font_size, labelpad=10, fontname="Times New Roman")
plt.ylabel('Normalized Intensity', fontsize=font_size, labelpad=10, fontname="Times New Roman")
plt.xticks([0,60,120,180,240], fontsize=font_size, fontname="Times New Roman")
plt.yticks(fontsize=font_size, fontname="Times New Roman")
#plt.savefig('C:/Users/Brendan/Desktop/%s_%i_timeseries.pdf' % (date, wavelength), format='pdf', bbox_inches='tight')

       
##########
###########
##########        
        

### Compute pixel's power spectrum + fit models to spectrum + plot    
    
cube_shape = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength))
spectra_array = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))
## load in array of segment-averaged pixel FFTs
SPECTRA = spectra_array


num_freq = 149

freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script

time_step = 24

sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

                         
f = freqs  # frequencies
s = spectra_array[y0-1][x0-1]  # fourier power

   
# assign equal weights to all parts of the curve
df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
df2 = np.zeros_like(f)
df2[0:len(df)] = df
df2[len(df2)-1] = df2[len(df2)-2]
ds = df2
                                       
### fit data to models using SciPy's Levenberg-Marquart method

try:
    # initial guesses for fitting parameters
    M1_low = [-0.002, 0.3, -0.01]
    M1_high = [0.002, 6., 0.01]
    nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays

except RuntimeError:
    print("Error M1 - curve_fit failed")

except ValueError:
    print("Error M1 - inf/NaN ")
  
A, n, C = nlfit_l  # unpack fitting parameters


## fit data to combined power law plus gaussian component model
      
try:                                 
    M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
    M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
    #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
    
    # change method to 'dogbox' and increase max number of function evaluations to 3000
    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.43], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
    
except RuntimeError:
    print("Error M2 - curve_fit failed")

except ValueError:
    print("Error M2 - inf/NaN")

A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
print nlfit_gp

# unpack uncertainties in fitting parameters from diagonal of covariance matrix
dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]

  
# create model functions from fitted parameters
m1_fit = PowerLaw(f, A, n, C)
amp_scale = PowerLaw(np.exp(fp2), A, n, C)  # to extract the gaussian-amplitude scaling factor
m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)

nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
m2_param = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22) 
uncertainties = dA22, dn22, dC22, dP22, dfp22, dfw22  # do we want to keep a global array of uncertainties?

residsM22 = (s - m2_fit2)
chisqrM22 = ((residsM22/ds)**2).sum()
redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6)  

m2P_fit = PowerLaw(f, A22, n22, C22)  # only need if plotting
m2G_fit = Gauss(f, P22, fp22, fw22)  # only need if plotting                 

residsM2 = (s - m2_fit)
chisqrM2 = ((residsM2/ds)**2).sum()
redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)

residsM1 = (s - m1_fit)
chisqrM1 =  ((residsM1/ds)**2).sum()
redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       

residsM22 = (s - m2_fit2)
chisqrM22 = ((residsM22/ds)**2).sum()
redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6) 
    
f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))

amp_scale2 = PowerLaw(np.exp(fp22), A22, n22, C22)  # to extract the gaussian-amplitude scaling factor

# generate p-value heatmap
df1, df2 = 3, 6
p_val = ff.sf(f_test2, df1, df2)

r_val = pearsonr(m2_fit2, s)


fig = plt.figure(figsize=(12,10))
ax = plt.gca()  # get current axis -- to set colorbar 
plt.title(r'(d) Power Spectrum', y = 1.01, fontsize=25)
plt.ylim((10**-4.7,10**0))
plt.xlim((10**-4.,10**-1.3))
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax.tick_params(axis='both', which='major', pad=10)
plt.loglog(f,s,'k')
plt.loglog(f, m2_fit2, 'purple', label='M2 - Combined', linewidth=1.3)
plt.loglog(f, m2G_fit, 'g--', label='M2 - Gaussian', linewidth=1.3)
plt.xlabel('Frequency [Hz]', fontsize=font_size, labelpad=10)
plt.ylabel('Power', fontsize=font_size, labelpad=10)
plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')

plt.text(0.00763, 10**-0.75, r'$\beta$ = {0:0.2f} [min]'.format((1./np.exp(m2_param[4]))/60.), fontsize=font_size)
plt.text(0.00786, 10**-1., r'$r$ = {0:0.3g}'.format(r_val[0]), fontsize=font_size)
legend = ax.legend(loc='lower left', prop={'size':font_size}, labelspacing=0.35)
for label in legend.get_lines():
    label.set_linewidth(3.0)  # the legend line width
#plt.savefig('C:/Users/Brendan/Desktop/%s_%i_spectrum_gaussian.pdf' % (date, wavelength), format='pdf', bbox_inches='tight')
    
    
##########
###########
##########        


### generate visual image

titles_vis = ['Average', 'Middle-File']
names_vis = ['average', 'mid']

vis = np.load('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength))

v_min = np.percentile(vis[0],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
v_max = np.percentile(vis[0],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)      

x_ticks = [10,210,410,610,810,1010,1210]
y_ticks = [10,210,410,610,810,1010,1210]  
x_ind = [-600,-400,-200,0,200,400,600]
y_ind = [600,400,200,0,-200,-400,-600]   

#y0 = vis.shape[1]-800
#x0 = 800  
y0 = vis.shape[1]-y0
x0 = x0

fig = plt.figure(figsize=(12,10))
ax1 = plt.gca()
plt.subplots_adjust(right=0.875)  #to substitute for colorbar space
plt.title(r'(a) Visual Average', y = 1.01, fontsize=font_size)  # no date / wavelength
im = ax1.imshow(np.flipud(vis[0]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
plt.scatter(x0, y0, s=100, c='red')
plt.scatter(x0, y0, s=50, c='red')
plt.scatter(x0, y0, s=75, c='white', marker='*')
plt.xlim(0, vis[0].shape[1])
plt.ylim(vis[0].shape[0], 0)
plt.xticks(x_ticks,x_ind,fontsize=font_size)
plt.yticks(y_ticks,y_ind,fontsize=font_size)
ax1.tick_params(axis='both', which='major', pad=10)
rect = patches.Rectangle((775,403), 50, 50, color='red', fill=False, linewidth=2)  # 1600
#rect = patches.Rectangle((776,398), 50, 50, color='red', fill=False, linewidth=2)  # 1700
ax1.add_patch(rect)

#plt.savefig('C:/Users/Brendan/Desktop/%s_%i_visual_point.pdf' % (date, wavelength), format='pdf', bbox_inches='tight')


"""
###################################################
#####################################################
###################################################
"""



### Make 4-panel figure for paper: visual average w/ selected pixel, time series, spectra + fit, and Gaussian location

plt.rcParams["font.family"] = "Times New Roman"
font_size = 12  # set the font size to be used for all text - titles, tick marks, text, labels


## visual image w/ selected pixel

fig = plt.figure(figsize=(12,10))
ax1 = plt.gca()
ax1 = plt.subplot2grid((11,66),(0, 0), colspan=29, rowspan=5)  #to substitute for colorbar space
plt.title(r'(a) Visual Average', y = 1.01, fontsize=font_size)  # no date / wavelength
im = ax1.imshow(np.flipud(vis[0]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
#plt.scatter(x0, y0, s=100, c='red')
#plt.scatter(x0, y0, s=50, c='red')
#plt.scatter(x0, y0, s=75, c='white', marker='*')
plt.xlim(0, vis[0].shape[1])
plt.ylim(vis[0].shape[0], 0)
plt.xticks(x_ticks,x_ind,fontsize=font_size)
plt.yticks(y_ticks,y_ind,fontsize=font_size)
ax1.tick_params(axis='both', which='major', pad=10)
rect = patches.Rectangle((771,402), 50, 50, color='red', fill=False, linewidth=2)  # 1600 (775,403 - 1600)
#rect = patches.Rectangle((776,398), 50, 50, color='red', fill=False, linewidth=2)  # 1700
ax1.add_patch(rect)


## time series of selected pixel

ax2 = plt.gca()
ax2 = plt.subplot2grid((11,66),(0, 36), colspan=30, rowspan=5)  #to substitute for colorbar space
plt.title(r'(b) Time Series', y = 1.01, fontsize=font_size)
ax2.tick_params(axis='both', which='major', pad=10)
plt.xlim(0,240)
plt.plot(TIME/60., timeseries, 'k')
plt.xlabel(r'Time [min]', fontsize=font_size, labelpad=10, fontname="Times New Roman")
plt.ylabel('Normalized Intensity', fontsize=font_size, labelpad=10, fontname="Times New Roman")
plt.xticks([0,60,120,180,240], fontsize=font_size, fontname="Times New Roman")
plt.yticks(fontsize=font_size, fontname="Times New Roman")


## compute Gaussian location w/ mask for region

h_map = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))

# generate p-value heatmap + masked Gaussian component heatmaps
df1, df2 = 3, 6  # degrees of freedom for model M1, M2
p_val = ff.sf(h_map[6], df1, df2)


mask_thresh = 0.005  # significance threshold - masked above this value
   
p_mask = np.copy(p_val)
loc_mask = np.copy(h_map[4])    

# mask the Gaussian component arrays with NaNs if above threshold 
p_mask[p_val > mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
loc_mask[p_val > mask_thresh] = np.NaN

# for creating sunspot umbra + PPV contour overlays from 1600
v1600 = np.load('%s/DATA/Output/%s/1600/visual.npy'% (directory, date))  
v1600 = v1600[:,1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
p1600 = np.load('%s/DATA/Output/%s/1600/param.npy'% (directory, date)) 

h_map = h_map[:,:p1600.shape[1],:p1600.shape[2]]

p1600_val = ff.sf(p1600[6], df1, df2)
v_mask = np.copy(v1600[0])
   
v_mask[p1600_val < mask_thresh] = 1.  # invert mask, set equal to 1. -- so can make contour

# determine percentage of region masked 
count = np.count_nonzero(np.isnan(p_mask))   
total_pix = p_val.shape[0]*p_val.shape[1]
mask_percent = ((np.float(count))/total_pix)*100
                
loc_mask = (1./np.exp(loc_mask))/60.  # convert Gaussian location to minutes

h_map[4] = (1./(np.exp(h_map[4])))/60.
#h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
#h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
h_min = 3.5
h_max = 5.5
cmap = cm.get_cmap('jet_r', 10)

h_range = np.abs(h_max-h_min)
h_step = h_range / 10.
c_ticks = np.zeros((11))
for h in range(11):
    c_ticks[h] = h_min + h_step*h 


## Gaussian location heatmap

ax3 = plt.gca()
ax3 = plt.subplot2grid((11,66),(6, 0), colspan=30, rowspan=5)  #to substitute for colorbar space
plt.title(r'(c) Gauss. Loc. $\beta$ [min]; $p$ < %0.3f | f$_{masked}$ = %0.1f%s' % (mask_thresh, mask_percent, '%'), y = 1.01, fontsize=font_size)
im = ax3.imshow(np.flipud(loc_mask), cmap = cmap, vmin=h_min, vmax=h_max)
plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax3)  # set colorbar to heatmap axis
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax, format='%0.1f')
cbar.ax.tick_params(labelsize=font_size, pad=5) 
cbar.set_ticks(c_ticks)


## power spectrum + fit of selected pixel

ax4 = plt.gca()
ax4 = plt.subplot2grid((11,66),(6, 36), colspan=30, rowspan=5)  #to substitute for colorbar space
plt.title(r'(d) Power Spectrum', y = 1.01, fontsize=font_size)
plt.ylim((10**-4.7,10**0))
plt.xlim((10**-4.,10**-1.3))
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax4.tick_params(axis='both', which='major', pad=5)
plt.loglog(f,s,'k')
plt.loglog(f, m2_fit2, 'purple', label='M2 - Combined', linewidth=1.3)
plt.loglog(f, m2G_fit, 'g--', label='M2 - Gaussian', linewidth=1.3)
plt.xlabel('Frequency [Hz]', fontsize=font_size, labelpad=5)
plt.ylabel('Power', fontsize=font_size, labelpad=5)
plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')

plt.text(0.00763, 10**-0.75, r'$\beta$ = {0:0.2f} [min]'.format((1./np.exp(m2_param[4]))/60.), fontsize=font_size)
plt.text(0.00786, 10**-1., r'$r$ = {0:0.3g}'.format(r_val[0]), fontsize=font_size)
legend = ax4.legend(loc='lower left', prop={'size':font_size}, labelspacing=0.35)
for label in legend.get_lines():
    label.set_linewidth(3.0)  # the legend line width
    
#plt.savefig('C:/Users/Brendan/Desktop/%s_%i_summary_revA.pdf' % (date, wavelength), format='pdf', bbox_inches='tight', dpi=300)
    
    

#######
#########
#######

## Plot combined histogram for 1600 + 1700

wavelength = 1600

# load parameter array and visual images from file tree structure 
heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
visual = np.load('%s/DATA/Output/%s/%i/visual.npy'% (directory, date, wavelength))  

visual = visual[:,1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
h_map = heatmaps

p1600 = np.load('%s/DATA/Output/%s/1600/param.npy'% (directory, date)) 

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels


h_map[4] = (1./(np.exp(h_map[4])))/60.
h_min = np.percentile(h_map[4],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
h_max = np.percentile(h_map[4],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
cmap = cm.get_cmap('jet_r', 10)


flat_param = np.reshape(h_map[4], (h_map[4].shape[0]*h_map[4].shape[1]))

# calculate some statistics
mean = np.mean(flat_param)
sigma = np.std(flat_param)   

fig = plt.figure(figsize=(12,10))
plt.title('2012/01/11', y = 1.01, fontsize=font_size)  # no date / wavelength
plt.xlabel(r'Gauss. Loc. $\beta$ [min]', fontsize=font_size, labelpad=10)
plt.ylabel('Bin Count', fontsize=font_size, labelpad=10)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)

y, x, _ = plt.hist(flat_param, bins=205, range=(3.5,5.),color='green',alpha=0.0)
n=y[1:-2]
bins=x[1:-2]
elem = np.argmax(n)
bin_max = bins[elem]
plt.ylim(0, y.max()*1.1)
y, x, _ = plt.hist(flat_param, bins=205, range=(3.5,5.),color='green',alpha=0.7,label=r'1600 $\AA$ | Mode = %0.3f [min]' % bin_max)
   
plt.vlines(bin_max, 0, y.max()*1.1, color='black', linestyle='dashed', linewidth=3.)   


    
wavelength = 1700

# load parameter array and visual images from file tree structure 
heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
visual = np.load('%s/DATA/Output/%s/%i/visual.npy'% (directory, date, wavelength))  

visual = visual[:,1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
h_map = heatmaps

p1600 = np.load('%s/DATA/Output/%s/1600/param.npy'% (directory, date)) 

h_map[4] = (1./(np.exp(h_map[4])))/60.
h_min = np.percentile(h_map[4],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
h_max = np.percentile(h_map[4],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
cmap = cm.get_cmap('jet_r', 10)


flat_param = np.reshape(h_map[4], (h_map[4].shape[0]*h_map[4].shape[1]))

# calculate some statistics
mean = np.mean(flat_param)
sigma = np.std(flat_param)   


y, x, _ = plt.hist(flat_param, bins=205, range=(3.5,5.),color='purple',alpha=0.0)
n=y[1:-2]
bins=x[1:-2]
elem = np.argmax(n)
bin_max = bins[elem]
plt.ylim(0, y.max()*1.1)
y, x, _ = plt.hist(flat_param, bins=205, range=(3.5,5.),color='purple',alpha=0.7,label=r'1700 $\AA$ | Mode = %0.3f [min]' % bin_max)
    
plt.vlines(bin_max, 0, y.max()*1.1, color='black', linestyle='dashed', linewidth=3.)   
#plt.hlines(y[nearest],fwhm_min,fwhm_min+fwhm, linestyle='dashed', linewidth=2., color='white')
legend = plt.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
for label in legend.get_lines():
    label.set_linewidth(2.0)  # the legend line width
plt.xlim(3.5,5.)

#plt.savefig('C:/Users/Brendan/Desktop/1600_1700_%s_combined_hist.pdf' % date, format='pdf', bbox_inches='tight')