# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 13:47:34 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm
from scipy import stats
import scipy.signal
import jdcal
from astropy.time import Time
import datetime
from numpy.random import choice
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from pylab import *
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import h5py
from scipy import fftpack  # doesnt work in module when called here???
from astropy.convolution import convolve, Box1DKernel
from numpy.random import randn
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from timeit import default_timer as timer
from scipy.stats import f
import matplotlib.patches as patches
from scipy.stats.stats import pearsonr

from scipy import fftpack    

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)

fmt = '%Y%m%d'

directory = 'S:'

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels

#bstrp_param_arr = np.load('C:/Users/Brendan/Desktop/Inbox/1600_1700_bootstrap_params/1600_1700_bootstrap_params_array.npy')
bstrp_param_arr = np.load('C:/Users/Brendan/Desktop/1600_1700_bootstrap_params_array_NovB.npy')

#cal_param_arr = np.load('C:/Users/Brendan/Desktop/UV_calibrated_params.npy')
#cal_param_arr = np.load('C:/Users/Brendan/Desktop/UV_calibrated_params_NovB.npy')
cal_param_arr = np.load('C:/Users/Brendan/Desktop/UV_calibrated_params_Nov_M2.npy')


hist_stddev = np.load('C:/Users/Brendan/Desktop/1600_1700_gauss_loc_std_dev.npy')
#uncertainty = np.load('C:/Users/Brendan/Desktop/1600_1700_gauss_loc_mode_uncertainty.npy')
uncertainty = np.load('C:/Users/Brendan/Desktop/1600_1700_gauss_loc_mode_uncertainty_NovB.npy')

dates = bstrp_param_arr[0][9]
jul_dates = bstrp_param_arr[0][8]

greg_dates = []

for q in range(len(dates)):
    greg_date_temp = '%s/%s/%s' % (str(dates[q])[0:4], str(dates[q])[4:6], str(dates[q])[6:8])
    greg_dates = np.append(greg_dates, greg_date_temp)
    

raw1600 = bstrp_param_arr[0]    
raw1700 = bstrp_param_arr[1]    
cal1600 = cal_param_arr[0]
cal1700 = cal_param_arr[1]



#"""
#sunspot_data = np.load('C:/Users/Brendan/Desktop/Inbox/SDO_timespan_sunspot_data.npy')
sunspot_data = np.load('C:/Users/Brendan/Desktop/Inbox/SDO_timespan_sunspot_data_8_17.npy')

sunspot_data[0] /= np.max(sunspot_data[0])
sunspot_data[1] /= np.max(sunspot_data[1])
sunspot_data[2] /= np.max(sunspot_data[2])

p = 4

if p == 1:
    param1600_rev = cal1600[p,:]
    param1700_rev = cal1700[p,:]
    param_name = 'Power Law Index'
    save_name = 'index'
elif p == 4:
    param1600 = raw1600[p]
    param1700 = raw1700[p]
    param1600_rev = (1./np.exp(cal1600[p,:]))/60.
    param1700_rev = (1./np.exp(cal1700[p,:]))/60.
    param_name = 'Gaussian Location'
    save_name = 'gauss_loc'

offset = param1700_rev - param1600_rev
#p1600rev_interp = np.interp(sunspot_data[3],j_dates,param1600_rev)
#p1700rev_interp = np.interp(sunspot_data[3],j_dates,param1700_rev)

#p1600_interp /= np.max(p1600_interp)
#p1700_interp /= np.max(p1700_interp)

#r1600 = pearsonr(p1600_interp, sunspot_data[1])[0]  # calculate r-value correlation coefficient
#r1600_rev = pearsonr(p1600rev_interp, sunspot_data[1])[0]  # calculate r-value correlation coefficient
#r1700 = pearsonr(p1700_interp, sunspot_data[1])[0]  # calculate r-value correlation coefficient

# calculate some statistics
mean1600 = np.mean(param1600_rev)
sigma1600 = np.std(param1600_rev)
mean1700 = np.mean(param1700_rev)
sigma1700 = np.std(param1700_rev)
mean_offset = np.mean(offset)
sigma_offset = np.std(offset)

fig = plt.figure(figsize=(17,10))

ymin = np.min([param1600_rev, param1700_rev])*0.99
#ymax = np.max([centers1600_bstrp, centers1700_bstrp])*1.035
#ymax = np.max([param1600_rev, param1700_rev])*1.015
ymax = np.max([param1600_rev, param1700_rev])*1.03
#ymin = np.min([param1600_rev, param1700_rev])*0.80 # for distribution std_dev
#ymax = np.max([param1600_rev, param1700_rev])*1.20
   
ax2 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
#ax2.set_title(r'%s | Radio Flux' % param_name, y=1.01, fontsize=font_size)
#ax2.set_title(r'%s' % param_name, y=1.01, fontsize=font_size)
#ax2.set_title(r'%s | Calibrated' % param_name, y=1.01, fontsize=font_size)
#ax2.set_title(r'%s | Raw' % param_name, y=1.01, fontsize=font_size)
#ax2.set_title(r'%s' % param_name, y=1.01, fontsize=font_size)
ax2.set_title(r'Dominant Peak Oscillation', y=1.01, fontsize=font_size)

#ax2.scatter(jul_dates, param1700_rev, color='purple', linewidth=3., label=r'1700 $\AA$')
ax2.scatter(jul_dates, param1700, color='white', edgecolors='purple', marker='^', linewidth=2., s=100, alpha=0.9, label=r'1700 $\AA$ (Raw)')
#ax2.scatter(jul_dates, param1700, color='white', marker='^', linewidth=1.)
ax2.scatter(jul_dates, param1700_rev, color='purple', marker='o', linewidth=3., label=r'1700 $\AA$ (Calibrated)')
#ax2.scatter(jul_dates, param1700_rev, color='purple', linewidth=3.)
#l_1700 = ax2.plot(jul_dates, param1700_rev, linestyle='dashed', linewidth=3., color='purple', label=r'1700 $\AA$ | Mean = %0.3f | $\sigma$ = %0.3f' % (mean1700, sigma1700))
ax2.scatter(jul_dates, param1600, color='white', edgecolors='green',  marker='^', linewidth=2., s=100, alpha=0.9, label=r'1600 $\AA$ (Raw)')
#ax2.scatter(jul_dates, param1600, color='white', marker='^',  linewidth=1.)
ax2.scatter(jul_dates, param1600_rev, color='green', marker='o',  linewidth=3., label=r'1600 $\AA$ (Calibrated)')
#ax2.scatter(jul_dates, param1600_rev, color='green', linewidth=3.)
#l_1600 = ax2.plot(jul_dates, param1600_rev, linestyle='dashed', linewidth=3., color='green', label=r'1600 $\AA$ | Mean = %0.3f | $\sigma$ = %0.3f' % (mean1600, sigma1600))

#ax2.errorbar(jul_dates, param1700_rev, yerr=hist_stddev[1], ecolor='purple', elinewidth=1.5)
#ax2.errorbar(jul_dates, param1600_rev, yerr=hist_stddev[0], ecolor='green', elinewidth=1.5)
ax2.errorbar(jul_dates, param1700_rev, yerr=uncertainty[1], fmt='none', ecolor='purple', elinewidth=1.)
ax2.errorbar(jul_dates, param1600_rev, yerr=uncertainty[0], fmt='none', ecolor='green', elinewidth=1.)

ax2.set_ylabel('%s [min]' % param_name, fontsize=font_size)
ax2.set_ylim(ymin,ymax)
ax2.set_xlabel('Date', fontsize=font_size)
plt.xticks(jul_dates, greg_dates, rotation=60, fontsize=15)
plt.yticks(fontsize=15)
ax2.set_xlim(jul_dates[0]-60,jul_dates[-1]+60)
ax3 = ax2.twinx()
ax2.minorticks_on()
ax2.yaxis.grid(True, which='both')
ax2.xaxis.grid(False)
#l_offset = ax3.plot(jul_dates, offset, color='red', linestyle='solid', linewidth=2., label=r' Offset  | Mean = %0.3f | $\sigma$ = %0.3f' % (mean_offset, sigma_offset))

#labels = l_1700 + l_1600 + l_offset
#labls = [l.get_label() for l in labels]

legend2 = ax2.legend(loc='upper right', prop={'size':18}, labelspacing=0.35)
#legend2 = ax2.legend(labels, labls, loc='upper right', prop={'size':18}, labelspacing=0.35)
for label in legend2.get_lines():
    label.set_linewidth(4.0)  # the legend line width
    
#ax3.plot(sunspot_data[3], sunspot_data[0], color='red', linestyle='solid', linewidth=2.)
#ax3.plot(sunspot_data[3], sunspot_data[1], color='red', linestyle='solid', linewidth=2.)
#ax3.plot(sunspot_data[3], sunspot_data[2], color='red', linestyle='solid', linewidth=2.)
#ax3.bar(sunspot_data[3], sunspot_data[1], color='red', width=30, alpha=0.5)
#ax3.set_ylabel(r'Offset : 1700 $\AA$ - 1600 $\AA$', labelpad = 7, fontsize=font_size)
#ax3.set_ylim(0.1,0.2)
#plt.yticks(fontsize=15)
    
#plt.savefig('C:/Users/Brendan/Desktop/1600_1700_gauss_loc_raw_cal_errorbars_revB.pdf', format='pdf', bbox_inches='tight')



"""
index1600_bstrp = bstrp_param_arr[0][1]
index1700_bstrp = bstrp_param_arr[1][1]

gauss1600_bstrp = bstrp_param_arr[0][4]
gauss1700_bstrp = bstrp_param_arr[1][4]

index_diff = index1700_bstrp - index1600_bstrp    
gauss_diff = gauss1700_bstrp - gauss1600_bstrp

index_diff_rev = cal1700[:,1] - cal1600[:,1]
gauss_diff_rev = (1./np.exp(cal1700[:,4]))/60. - (1./np.exp(cal1600[:,4]))/60.

fig = plt.figure(figsize=(17,10))
#ymin = np.min([centers1600_bstrp, centers1700_bstrp])*0.98
#ymax = np.max([centers1600_bstrp, centers1700_bstrp])*1.02
ax2 = plt.gca()
#ax2 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
ax2.set_title(r'Difference in Power Law Index & Difference in Gaussian Location', y=1.01, fontsize=font_size)

ax2.scatter(jul_dates,gauss_diff, color='red', linewidth=3.)
ax2.plot(jul_dates,gauss_diff,linestyle='solid', linewidth=3., color='red', label=r'Gauss. Loc. (Original)')
ax2.scatter(j_dates,gauss_diff_rev, color='red', linewidth=3.)
ax2.plot(j_dates,gauss_diff_rev,linestyle='dashed', linewidth=3., color='red', label=r'Gauss. Loc. (Revised)')

ax2.scatter(jul_dates,index_diff, color='black', linewidth=3.)
ax2.plot(jul_dates,index_diff,linestyle='solid', linewidth=3., color='black', label=r'Index (Original)')
ax2.scatter(j_dates,index_diff_rev, color='black', linewidth=3.)
ax2.plot(j_dates,index_diff_rev,linestyle='dashed', linewidth=3., color='black', label=r'Index (Revised)')

ax2.set_ylabel(r'1700 $\AA$ - 1600 $\AA$', fontsize=font_size)
ax2.set_xlabel('Date', fontsize=font_size)
plt.xticks(jul_dates, greg_dates, rotation=60)
ax2.set_xlim(jul_dates[0]-26,jul_dates[-1]+7)
ax2.set_ylim(0.02,0.21)
#ax2.hlines(0.15,jul_dates[0]-26,jul_dates[-1]+7, linestyle='dashed', linewidth=2., color='blue')
#ax2.hlines(0.075,jul_dates[0]-26,jul_dates[-1]+7, linestyle='dashed', linewidth=2., color='blue')
ax2.minorticks_on()
ax2.yaxis.grid(True, which='both')
ax2.xaxis.grid(False)
#ax3 = ax2.twinx()
legend2 = ax2.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
for label in legend2.get_lines():
    label.set_linewidth(4.0)  # the legend line width

#plt.savefig('C:/Users/Brendan/Desktop/1600_1700_param_diff_calibrated_grid_rev.pdf', format='pdf', bbox_inches='tight')
"""