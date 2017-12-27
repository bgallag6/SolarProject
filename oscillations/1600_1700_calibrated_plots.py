# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:47:41 2017

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

#solar_activity_titles = ['IMF Magnitude Avg, nT', 'Bz, GSM, nT', 'Bz, GSE, nT', 'Proton Temperature, K', 'Proton Density, n/cc', 'Alpha/Proton Density Ratio', 'Flow Pressure, nPa', 'Ey - Electric Field, mV/m', 'Plasma Beta', 'Alfven Mach Number', 'Magnetosonic Mach Number', 'R Sunspot Number (new version)', 'Dst Index, nT', 'Solar index F10.7 ( sfu = 10-22.m-2.Hz-1)', 'AE Index, nT', 'AL Index, nT', 'Lyman Alpha Solar Index (1e11 photons/cm^2/sec)']
solar_activity_fnames = ['imf_Bfield', 'Bz_GSM', 'Bz_GSE', 'sw_proton_temp', 'sw_proton_density', 'alpha_proton_ratio', 'flow_pressure', 'Efield', 'plasma_beta', 'alfen_mach_number', 'magnetosonic_mach_number', 'sunspot_number', 'dst_index', 'f10_7_index', 'AE_index', 'AL_index', 'Lyman_alpha']
solar_activity_titles = ['Bz, GSE, nT', 'Bz, GSM, nT', 'Plasma Beta', 'R Sunspot Number (new version)', 'Solar index F10.7 ( sfu = 10-22.m-2.Hz-1)', 'Lyman Alpha Solar Index (1e11 photons/cm^2/sec)']

bstrp_param_arr = np.load('C:/Users/Brendan/Desktop/1600_1700_bootstrap_params_array_NovB.npy')

dates = bstrp_param_arr[0][9]
jul_dates = bstrp_param_arr[0][8]

greg_dates = []

for q in range(len(dates)):
    greg_date_temp = '%s/%s/%s' % (str(dates[q])[0:4], str(dates[q])[4:6], str(dates[q])[6:8])
    greg_dates = np.append(greg_dates, greg_date_temp)
    
    
time_step = 24  # 12-second cadence for the others
n_segments = 2 
n_freqs = 149
freq_size = ((n_freqs)*2) + 1  # determined from FFT-averaging script
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

uncertainty = np.load('C:/Users/Brendan/Desktop/1600_1700_gauss_loc_mode_uncertainty_NovB.npy')

#cal1600 = np.load('C:/Users/Brendan/Desktop/UV_calibrated_params_NovC.npy')[0]
#cal1700 = np.load('C:/Users/Brendan/Desktop/UV_calibrated_params_NovC.npy')[1]
cal1600 = np.load('C:/Users/Brendan/Desktop/UV_calibrated_params_Nov_M2.npy')[0]
cal1700 = np.load('C:/Users/Brendan/Desktop/UV_calibrated_params_Nov_M2.npy')[1]

#sunspot_data = np.load('C:/Users/Brendan/Desktop/Inbox/SDO_timespan_sunspot_data.npy')
#sunspot_data = np.load('C:/Users/Brendan/Desktop/Inbox/SDO_timespan_sunspot_data_8_17.npy')
#sunspot_data = np.load('C:/Users/Brendan/Desktop/Inbox/SDO_timespan_sunspot_data.npy')
#sunspot_data = np.load('C:/Users/Brendan/Desktop/solar_activity_full_curves.npy')
sunspot_data = np.load('C:/Users/Brendan/Desktop/solar_activity_full_curves_Nov.npy')

    
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

    
for i in range(5,6):
    
    sunspot_compare = i
    days_smooth = 27
    days_trim = days_smooth/2

    sunspot_data[sunspot_compare] = smooth(sunspot_data[sunspot_compare],days_smooth)
    s_data = sunspot_data[sunspot_compare][days_smooth:-(days_trim+1)]  # trim off smoothing edge effects (27 days from start because that is first point)
    #s_data /= np.max(s_data)  # normalize
    s_date = sunspot_data[6][days_smooth:-(days_trim+1)]

    
    p = 4
       
    centers1600_bstrp = bstrp_param_arr[0][p]
    centers1700_bstrp = bstrp_param_arr[1][p]
    
    if p == 1:
        param1600_rev = cal1600[p]
        param1700_rev = cal1700[p]
        param_name = 'Power Law Index'
        save_name = 'index'
    elif p == 3:
        param1600_rev = cal1600[p]
        param1700_rev = cal1700[p]
        param_name = 'Gaussian Amplitude'
        save_name = 'gauss_amp'
    elif p == 4:
        param1600_rev = (1./np.exp(cal1600[p]))/60.
        param1700_rev = (1./np.exp(cal1700[p]))/60.
        param_name = 'Gaussian Location'
        save_name = 'gauss_loc'
    elif p == 5:
        param1600_rev = cal1600[p]
        param1700_rev = cal1700[p]
        param_name = 'Gaussian Width'
        save_name = 'gauss_wid'
    
    #"""
    ### interpolate params onto solar activty indicators
    p1600_interp = np.interp(s_date,jul_dates,centers1600_bstrp)
    p1700_interp = np.interp(s_date,jul_dates,centers1700_bstrp)
    p1600rev_interp = np.interp(s_date,jul_dates,param1600_rev)
    p1700rev_interp = np.interp(s_date,jul_dates,param1700_rev)
    
    ## calculate r-value correlation coefficients
    r1600 = pearsonr(p1600_interp, s_data)[0]  
    r1700 = pearsonr(p1700_interp, s_data)[0]
    r1600_rev = pearsonr(p1600rev_interp, s_data)[0]
    r1700_rev = pearsonr(p1700rev_interp, s_data)[0]
    
    #p1600_interp /= np.max(p1600_interp)
    #p1700_interp /= np.max(p1700_interp)
    #"""
    
    """
    ### for interpolating solar activty indicator onto params
    sdata_interp = np.interp(jul_dates,s_date,s_data)
    
    r1600_rev = pearsonr(param1600_rev, sdata_interp)[0]
    r1700_rev = pearsonr(param1700_rev, sdata_interp)[0]
    """
    
    fig = plt.figure(figsize=(17,10))
    
    ymin = np.min([centers1600_bstrp, centers1700_bstrp])*0.98
    #ymax = np.max([centers1600_bstrp, centers1700_bstrp])*1.035
    ymax = np.max([centers1600_bstrp, centers1700_bstrp])*1.025
       
    ax2 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
    #ax2.set_title(r'%s | Radio Flux' % param_name, y=1.01, fontsize=font_size)
    #ax2.set_title(r'%s' % param_name, y=1.01, fontsize=font_size)
    #ax2.set_title(r'%s [%i-day Smooth]' % (solar_activity_titles[sunspot_compare], days_smooth), y=1.01, fontsize=font_size)
    #ax2.set_title(r'Gaussian Location', y=1.01, fontsize=font_size)
    #ax2.set_title(r'Gaussian Location | %s [%i-day Smooth]' % (solar_activity_titles[sunspot_compare], days_smooth), y=1.01, fontsize=font_size)
    
    ax2.errorbar(jul_dates, param1700_rev, yerr=uncertainty[1], ecolor='purple', elinewidth=1.5)
    ax2.errorbar(jul_dates, param1600_rev, yerr=uncertainty[0], ecolor='green', elinewidth=1.5)    
    
    #ax2.scatter(jul_dates,centers1600_bstrp, color='green', linewidth=3.)
    #ax2.scatter(jul_dates,centers1700_bstrp, color='purple', linewidth=3.)
    #ax2.plot(jul_dates,centers1700_bstrp,linestyle='dashed', linewidth=3., color='purple', label=r'1700 $\AA$ (Raw)')
    #ax2.plot(j_dates, param1700_rev, linestyle='solid', linewidth=3., color='purple', label=r'1700 $\AA$ (Calibrated)')
    ax2.plot(jul_dates, param1600_rev, linestyle='solid', linewidth=3., color='green', label=r'1600 $\AA$ | $r$-Value = %0.2f' % r1600_rev)
    #ax2.plot(jul_dates,centers1600_bstrp,linestyle='dashed', linewidth=3., color='green', label=r'1600 $\AA$ (Raw)')
    #ax2.plot(j_dates, param1600_rev, linestyle='solid', linewidth=3., color='green', label=r'1600 $\AA$ (Calibrated)')
    ax2.scatter(jul_dates, param1700_rev, color='purple', linewidth=3.)
    ax2.scatter(jul_dates, param1600_rev, color='green', linewidth=3.)
    ax2.plot(jul_dates, param1700_rev, linestyle='solid', linewidth=3., color='purple', label=r'1700 $\AA$ | $r$-Value = %0.2f' % r1700_rev)
    
    ax2.set_ylabel('%s [min]' % param_name, fontsize=font_size)
    ax2.set_ylim(ymin,ymax)
    ax2.set_xlabel('Date', fontsize=font_size)
    plt.xticks(jul_dates, greg_dates, rotation=60, fontsize=15)
    ax2.set_xlim(jul_dates[0]-60,jul_dates[-1]+60)
    plt.yticks(fontsize=15)
    ax3 = ax2.twinx()
    ax2.minorticks_on()
    ax2.yaxis.grid(True, which='both')
    ax2.xaxis.grid(False)
       
    legend2 = ax2.legend(loc='upper right', prop={'size':18}, labelspacing=0.35)
    for label in legend2.get_lines():
        label.set_linewidth(4.0)  # the legend line width
           
    ax3.plot(s_date, s_data, color='red', linestyle='solid', linewidth=2., alpha=0.5)
    #ax3.set_ylabel(r'Lyman Alpha [27-d Smooth, Normalized]', labelpad = 7, fontsize=font_size)
    #ax3.set_ylabel(r'Lyman Alpha [10$^{11}$ photons/cm$^2$/sec] | 27-d Smooth', labelpad = 7, fontsize=font_size)
    ax3.set_ylabel(r'Ly-$\alpha$ [10$^{11}$ photons/cm$^2$/sec] | 27-d Smooth', labelpad = 7, fontsize=font_size)
    #ax3.plot(j_dates, sdata_interp, color='red', linestyle='solid', linewidth=2., alpha=0.55)
    plt.yticks(fontsize=15)
    
    """
    ### no smoothing
    if i == 1: 
        ax3.set_ylim(-0.7,0.7)
    elif i == 4:
        ax3.set_ylim(0., 0.035)
    elif i == 5:
        ax3.set_ylim(0., 0.02)
    elif i == 6:
        ax3.set_ylim(0., 0.15)
    elif i == 7:
        ax3.set_ylim(-0.6, 0.6)
    elif i == 8:
        ax3.set_ylim(0., 0.035)
    elif i == 9:
        ax3.set_ylim(0., 0.035)
    elif i == 10:
        ax3.set_ylim(0.01, 0.12)
    elif i == 12:
        ax3.set_ylim(-5., 2.)
    elif i == 13:
        ax3.set_ylim(0.06, 0.165)
    elif i == 14:
        ax3.set_ylim(0., 0.1)
    elif i == 15:
        ax3.set_ylim(-0.006, 0.003)
    """
    
    #"""
    ### 27-day smoothing    
    if i == 0: 
        ax3.set_ylim(0.04,0.12)
    elif i == 1: 
        ax3.set_ylim(-0.025,0.02)
    elif i == 2: 
        ax3.set_ylim(0.015,0.045)
    elif i == 3:
        ax3.set_ylim(0., 1.1)
    elif i == 4:
        ax3.set_ylim(0.3, 1.1)
    #elif i == 5:
        #ax3.set_ylim(0.7, 1.05)
    #"""
    
    #plt.savefig('C:/Users/Brendan/Desktop/Inbox/1600_1700_bootstrap_params/1600_1700_bootstrap_%s_full_26_rval.pdf' % (names[c]), format='pdf', bbox_inches='tight')
    #plt.savefig('C:/Users/Brendan/Desktop/1600_1700_gauss_loc_lyman_alpha_%ismooth_revB.pdf' % (days_smooth), format='pdf', bbox_inches='tight')
    #plt.savefig('C:/Users/Brendan/Desktop/1600_1700_%s_%ismooth_%i_revA.pdf' % (save_name, days_smooth,i), format='pdf', bbox_inches='tight')
    #plt.savefig('C:/Users/Brendan/Desktop/1600_1700_gauss_loc_raw_calibrated.pdf', format='pdf', bbox_inches='tight')
    #plt.close()

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