# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 19:53:04 2017

@author: Brendan
"""

import numpy as np
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
from scipy import stats
from numpy.random import choice
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
from scipy import stats
from numpy.random import choice

n_samples = 100000
n_loops = 2000
#n_samples = 10000
#n_loops = 1000

fmt = '%Y%m%d'

directory = 'S:'
#dates = ['20101028','20110207', '20110601', '20110909', '20111030', '20120111', '20120319', '20120621', '20121227', '20130618', '20131015', '20140130', '20140601', '20140724', '20141205', '20150218', '20150530', '20150915', '20151231', '20160319', '20160602', '20160917', '20170111', '20170603', '20170803', '20170830']
dates = ['20171104']

titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'Power Law Index $n$', r'Power Law Tail - $C$', r'Gaussian Amplitude [flux] - $\alpha$', r'Gauss. Loc. $\beta$ [min]', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', r'$r$-Value: Correlation Coefficient', r'Rollover Period $T_r$ [min]'] # 10-param-list
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'r_value', 'roll_freq']

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels
    
greg_date = []
jul_date = []

for q in range(len(dates)):
    greg_date_temp = '%s/%s/%s' % (dates[q][0:4], dates[q][4:6], dates[q][6:8])
    greg_date = np.append(greg_date, greg_date_temp)
    dt = datetime.datetime.strptime(dates[q][0:10],fmt)
    jul_date_temp = sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) + dt.hour/24.
    jul_date = np.append(jul_date, jul_date_temp)
   
fig_width = 10+2  # works better for 20130626 (with no x/y labels)
fig_height = 10  # works better for 20130626

sigmas1600_bstrp = np.zeros((len(dates)))
centers1600_bstrp = np.zeros((len(dates)))
sigmas1700_bstrp = np.zeros((len(dates)))
centers1700_bstrp = np.zeros((len(dates)))

bstrp_param_arr = np.zeros((2,10,len(dates)))  # wavelength, parameter+date, date

bstrp_param_arr[0][8] = jul_date
bstrp_param_arr[0][9] = dates
bstrp_param_arr[1][8] = jul_date
bstrp_param_arr[1][9] = dates

for c in range(8):
#for c in range(1):  # params
    
    for i in range(len(dates)):
    #for i in range(3):  # dates
    
        fig = plt.figure(figsize=(20,10))
       
        ax1 = plt.subplot2grid((1,2),(0, 0), colspan=1, rowspan=1)
        ax2 = plt.subplot2grid((1,2),(0, 1), colspan=1, rowspan=1)
        plt.suptitle(r'%s' % titles[c], y=0.99, fontsize=font_size)
        ax1.set_title(r'1600 $\AA$', y=1.01, fontsize=font_size)
        ax2.set_title(r'1700 $\AA$', y=1.01, fontsize=font_size)
        
        date = dates[i]
        
        wavelength = 1600  
    
        # load parameter array and visual images from file tree structure 
        heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
        
        heatmaps[4] = (1./(np.exp(heatmaps[4])))/60.
        
        if i == 2 or i == 9 or i == 14:
            if c == 0:
                h_min = np.percentile(heatmaps[c],5)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
                h_max = np.percentile(heatmaps[c],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        else:
            h_min = np.percentile(heatmaps[c],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(heatmaps[c],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            
        # Make data 1-D
        param_1d = np.reshape(heatmaps[c], (heatmaps[c].shape[0]*heatmaps[c].shape[1]))
        
        all_max_vals=np.empty(n_loops) # Hold the results
        
        # LOOP SOME LARGE NUMBER OF TIMES
        for k in range(n_loops):
             sample = choice(param_1d, n_samples, replace=True) # GET RANDOM SAMPLE OF LOC_MASK AND FIND MODE
             #y,x=np.histogram(sample,bins=100,range=(h_min, h_max))
             y,x=np.histogram(sample,bins=200,range=(h_min, h_max))
             n=y[1:-2]
             bins=x[1:-2]
             elem = np.argmax(n)
             all_max_vals[k] = bins[elem]
             
        sigmas1600_bstrp[i] = np.std(all_max_vals)
        centers1600_bstrp[i] = np.mean(all_max_vals)   
        
        y, x, _ = ax1.hist(param_1d, range=(h_min, h_max), bins=200)
        n=y[1:-2]
        bins=x[1:-2]
        elem = np.argmax(n)
        mode = bins[elem]
        ax1.vlines(mode, 0,y.max()*1., color='red', linestyle='dashed', linewidth=2., label='Mode = %0.6f' % mode) 
        ax1.legend()
    
    
        wavelength = 1700
    
        # load parameter array and visual images from file tree structure 
        heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))    
        
        heatmaps[4] = (1./(np.exp(heatmaps[4])))/60.
        
        h_min = np.percentile(heatmaps[c],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(heatmaps[c],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            
        # Make data 1-D
        param_1d = np.reshape(heatmaps[c], (heatmaps[c].shape[0]*heatmaps[c].shape[1]))
     
        all_max_vals=np.empty(n_loops) # Hold the results
        
        # LOOP SOME LARGE NUMBER OF TIMES
        for k in range(n_loops):
             sample = choice(param_1d, n_samples, replace=True) # GET RANDOM SAMPLE OF LOC_MASK AND FIND MODE
             #y,x=np.histogram(sample,bins=100,range=(h_min, h_max))
             y,x=np.histogram(sample,bins=200,range=(h_min, h_max))
             n=y[1:-2]
             bins=x[1:-2]
             elem = np.argmax(n)
             all_max_vals[k] = bins[elem]
             
        sigmas1700_bstrp[i] = np.std(all_max_vals)
        centers1700_bstrp[i] = np.mean(all_max_vals)
        
        y, x, _ = ax2.hist(param_1d, range=(h_min,h_max), bins=200)
        n=y[1:-2]
        bins=x[1:-2]
        elem = np.argmax(n)
        mode = bins[elem]
        ax2.vlines(mode, 0,y.max()*1., color='red', linestyle='dashed', linewidth=2., label='Mode = %0.6f' % mode)
        ax2.legend()
        
        #plt.savefig('C:/Users/Brendan/Desktop/Inbox/1600_1700_bootstrap_params/1600_1700_bootstrap_%s_%s.pdf' % (names[c], date), format='pdf', bbox_inches='tight')
        plt.close()
    
    
    bstrp_param_arr[0][c] = centers1600_bstrp
    bstrp_param_arr[1][c] = centers1700_bstrp
    
    
    fig = plt.figure(figsize=(17,10))
    
    ymin = np.min([centers1600_bstrp, centers1700_bstrp])*0.98
    ymax = np.max([centers1600_bstrp, centers1700_bstrp])*1.02
       
    ax2 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
    ax2.set_title(r'1600 $\AA$ vs 1700 $\AA$: %s' % titles[c], y=1.01, fontsize=font_size)
    #ax2.scatter(jul_date,sigmas1600, color='green', linewidth=2.)
    #ax2.scatter(jul_date,sigmas1700, color='purple', linewidth=2.)
    ax2.scatter(jul_date,centers1600_bstrp, color='green', linewidth=2.)
    ax2.scatter(jul_date,centers1700_bstrp, color='purple', linewidth=2.)
    ax2.plot(jul_date,centers1600_bstrp,linestyle='dashed', linewidth=2., color='green', label=r'1600 $\AA$')
    ax2.plot(jul_date,centers1700_bstrp,linestyle='dashed', linewidth=2., color='purple', label=r'1700 $\AA$')
    ax2.set_ylabel('%s' % titles[c], fontsize=font_size)
    ax2.set_ylim(ymin,ymax)
    ax2.set_xlabel('Date', fontsize=font_size)
    plt.xticks(jul_date, greg_date, rotation=60)
    ax2.set_xlim(jul_date[0]-60,jul_date[-1]+60)
    #ax3 = ax2.twinx()
    #ax3.scatter(jul_date,sigmas1600, color='green', linewidth=2.)
    #ax3.scatter(jul_date,sigmas1700, color='purple', linewidth=2.)
    #ax3.plot(jul_date,sigmas1600,linestyle='dashed', linewidth=2., color='green', label='1600: Std. Dev.')
    #ax3.plot(jul_date,sigmas1700,linestyle='dashed', linewidth=2., color='purple', label='1700: Std. Dev')
    #ax3.set_ylabel('Sigma', fontsize=font_size)
    legend2 = ax2.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
    for label in legend2.get_lines():
        label.set_linewidth(2.0)  # the legend line width
        
    #plt.savefig('C:/Users/Brendan/Desktop/Inbox/1600_1700_bootstrap_params/1600_1700_bootstrap_%s_full.pdf' % (names[c]), format='pdf', bbox_inches='tight')
    plt.close()
    
#np.save('C:/Users/Brendan/Desktop/Inbox/1600_1700_bootstrap_params/1600_1700_bootstrap_params_array.npy', bstrp_param_arr)

"""
jul_dates_add = []

fmt = '%Y%m%d'

for r in range(len(dates_add)/2):
    dt = datetime.datetime.strptime(dates_add[(2*r)][0:10],fmt)
    jul_date_add = sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) + dt.hour/24.
    params_add[(2*r)] = np.append(params_add[(2*r)], jul_date_add)
    params_add[(2*r)] = np.append(params_add[(2*r)], float(dates_add[(2*r)]))
    params_add[(2*r)+1] = np.append(params_add[(2*r)+1], jul_date_add)
    params_add[(2*r)+1] = np.append(params_add[(2*r)+1], float(dates_add[(2*r)+1]))
    param_double = np.vstack((params_add[(2*r)], params_add[(2*r)+1]))
    bstrp_param_arr = np.insert(bstrp_param_arr, dates_ind[r], param_double, axis=2)  # inserts before designated index
    
    
#np.save('C:/Users/Brendan/Desktop/Inbox/1600_1700_bootstrap_params/1600_1700_bootstrap_params_array_26.npy', bstrp_param_arr)


#super_arr = np.zeros((2,10,len(dates)+3))
#super_arr[:,:,:len(dates)] = bstrp_param_arr
""" 
""" 
add_param_arr = np.zeros((2,10))

params_add = [0.000015, 1.022069, -0.000991, 0.019379, 4.056588, 0.239624, 83.153908, 4.921555]
date_add = '20140724'
dt = datetime.datetime.strptime(date_add[0:10],fmt)
jul_date_add = sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) + dt.hour/24.
ax2.scatter(jul_date_add, params_add[c], color='green')

add_param_arr[0][0:8] = params_add
add_param_arr[0][8] = jul_date_add
add_param_arr[0][9] = float(date_add)

super_arr[:,:,len(dates)] = add_param_arr

add_param_arr = np.zeros((2,10))

params_add = [0.000014, 1.027291, -0.000966, 0.019341, 4.079610, 0.236380, 83.530254, 4.863817]
date_add = '20150915'
dt = datetime.datetime.strptime(date_add[0:10],fmt)
jul_date_add = sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) + dt.hour/24.
ax2.scatter(jul_date_add, params_add[c], color='green')

add_param_arr[0][0:8] = params_add
add_param_arr[0][8] = jul_date_add
add_param_arr[0][9] = float(date_add)

super_arr[:,:,len(dates)+1] = add_param_arr

add_param_arr = np.zeros((2,10))

params_add = [0.000012, 1.068064, -0.001046, 0.018614, 4.173472, 0.230649, 73.292555, 4.571051]
date_add = '20150915'
dt = datetime.datetime.strptime(date_add[0:10],fmt)
jul_date_add = sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) + dt.hour/24.
ax2.scatter(jul_date_add, params_add[c], color='purple')

add_param_arr[1][0:8] = params_add
add_param_arr[1][8] = jul_date_add
add_param_arr[1][9] = float(date_add)

super_arr[:,:,len(dates)+2] = add_param_arr
"""