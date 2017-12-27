# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 10:40:16 2017

@author: Brendan
"""

from scipy.stats.stats import pearsonr 
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import scipy.signal
import jdcal
from astropy.time import Time
import datetime
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from pylab import *
from scipy.interpolate import interp1d
import scipy.misc
from scipy import fftpack  # doesnt work in module when called here???
  

# define Power-Law-fitting function (Model M1)
def PowerLaw(f2, A, n, C):
    return A*f2**-n + C
    
# define Gaussian-fitting function
def Gauss(f2, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f2))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)

fmt = '%Y%m%d'

directory = 'S:'

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels

#bstrp_param_arr = np.load('C:/Users/Brendan/Desktop/Inbox/1600_1700_bootstrap_params/1600_1700_bootstrap_params_array.npy')
bstrp_param_arr = np.load('C:/Users/Brendan/Desktop/1600_1700_bootstrap_params_array_NovB.npy')

dates = bstrp_param_arr[0][9]
jul_dates = bstrp_param_arr[0][8]

greg_dates = []

for q in range(len(dates)):
    greg_date_temp = '%s/%s/%s' % (str(dates[q])[0:4], str(dates[q])[4:6], str(dates[q])[6:8])
    greg_dates = np.append(greg_dates, greg_date_temp)
    
    
## determine frequencies to evaluate    
time_step = 24  # cadence for 1600/1700 A
n_freqs = 149
freq_size = ((n_freqs)*2) + 1
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

resids1600 = np.zeros((149))
resids1700 = np.zeros((149))
UV_slopes = np.zeros((2,149))
UV_intercept = np.zeros((2,149))  # [0] = 1600, [1] = 1700

## loop through wavelengths 
for w in range(2):
    
    ## loop through dates
    for i in range(bstrp_param_arr.shape[2]):
    
        A, n, C, P, fp, fw = bstrp_param_arr[w,:6,i]  # extract averaged parameters for each dataset
    
        fp = np.log((1./(fp*60.)))  # change Gaussian location to frequency from period
        
        #print A, n, C, P, fp, fw
        
        m2_spec = GaussPowerBase(freqs, A,n,C,P,fp,fw)  # generate averaged M2 spectra
        #m1_spec = PowerLaw(freqs, A, n, C)  # generate averaged M1 spectra
        #m2P_fit = PowerLaw(freqs, A, n, C)
        #m2G_fit = Gauss(freqs, P, fp, fw)
        #m2P = PowerLaw(freqs, bstrp_param_arr[w,0,0], bstrp_param_arr[w,1,0],bstrp_param_arr[w,2,0])  # generate baseline spectra (from initial dataset)
        m2P = GaussPowerBase(freqs, bstrp_param_arr[w,0,0], bstrp_param_arr[w,1,0],bstrp_param_arr[w,2,0], bstrp_param_arr[w,3,0], bstrp_param_arr[w,4,0],bstrp_param_arr[w,5,0])  # generate baseline spectra (from initial dataset)
        
        resids_temp = m2_spec - m2P  # calculate residuals between given M2 spectra and baseline
        #resids_temp = m1_spec - m2P  # calculate residuals between given M1 spectra and baseline
        
        ## add residuals for each dataset to array
        if w == 0:
            resids1600 = np.vstack((resids1600, resids_temp))
        elif w == 1:
            resids1700 = np.vstack((resids1700, resids_temp))
        
    
    ## for each frequency, extract residual from each date + calculate slope of best-fit line
    if w == 0:    
        resids1600 = resids1600[1:,:]   
        #plt.figure()
        for k in range(len(freqs)):
            f_resids = resids1600[:,k]
            m1, b1 = np.polyfit(jul_dates, f_resids, 1) 
            #print k, m1, b1
            UV_slopes[0,k] = m1
            UV_intercept[0,k] = b1
            #plt.plot(jul_dates, f_resids)
            #plt.plot(jul_dates, jul_dates*m1 + b1)
          
    elif w == 1:    
        resids1700 = resids1700[1:,:]   
        #plt.figure()
        for k in range(len(freqs)):
            f_resids = resids1700[:,k]
            m1, b1 = np.polyfit(jul_dates, f_resids, 1) 
            #print k, m1, b1
            UV_slopes[1,k] = m1
            UV_intercept[1,k] = b1
            #plt.plot(jul_dates, f_resids)
            #plt.plot(jul_dates, jul_dates*m1 + b1)

#np.save('C:/Users/Brendan/Desktop/UV_calibration_slopes.npy', UV_slopes)      


## to look at specific frequency interval    
#f_start = int(0)
#f_stop = int(149)  


"""
### Plotting

## plot residual curves by frequency
fig = plt.figure(figsize=(12,8))
plt.title(r'1600 $\AA$ | Slopes of Residuals by Frequency', fontsize=23, y=1.01)
for k in range(len(freqs)):
    plot_freq = freqs[k]*1000.
    plt.plot(jul_dates, resids1600[:,k], label='%0.2f mHz' % plot_freq)
    #plt.legend(loc='upper left', prop={'size':20}, labelspacing=0.35)
#plt.savefig('C:/Users/Brendan/Desktop/1600_spectra_resids_by_freq.pdf', format='pdf', bbox_inches='tight')


## plot slopes of residual curves by frequency
m2, b2 = np.polyfit(freqs, UV_slopes[0], 1) 
fig = plt.figure(figsize=(12,8))
plt.title(r'1600 $\AA$ | Slopes of Residuals by Frequency', fontsize=23, y=1.01)
plt.plot(freqs, UV_slopes[0]) 
#plt.plot(freqs[f_start:f_stop], UV_slopes[0,k][f_start:f_stop]) 
#plt.plot(freqs[f_start:f_stop], freqs[f_start:f_stop]*m2[f_start:f_stop]+b2[f_start:f_stop])        
#plt.loglog(freqs[f_start:f_stop], UV_slopes[0,k][f_start:f_stop]) 
#plt.loglog(freqs[f_start:f_stop], freqs[f_start:f_stop]*m2[f_start:f_stop]+b2[f_start:f_stop])      
#plt.savefig('C:/Users/Brendan/Desktop/1600_spectra_resids_slopes_by_freq_normal_rev.pdf', format='pdf', bbox_inches='tight')


f_stop = 149

## plot power-law portion of slopes of residual curves by frequency
nlfit_gp1, nlpcov_gp1 = scipy.optimize.curve_fit(PowerLaw, freqs[1:f_stop], UV_slopes[0][1:f_stop])
a1, b1, c1 = nlfit_gp1  # unpack fitting parameters     
fit1 = PowerLaw(freqs[1:f_stop], a1, b1, c1)

fig = plt.figure(figsize=(12,8))
plt.title(r'1600 $\AA$ | Slopes of Residuals by Frequency', fontsize=23, y=1.01)
plt.loglog(freqs[1:f_stop], UV_slopes[0][1:f_stop], linewidth=1.)  
plt.loglog(freqs[1:f_stop], fit1, linewidth=2., linestyle='dashed') 
plt.ylim(10**-10,10**-5)
plt.xlim(10**-4, 10**-1) 
plt.text(10**-2.5,10**-6.5,r'f($\nu$) = %0.2e$\cdot$$\nu$$^{-%0.3f}$' % (a1,b1), fontsize=23)  
#plt.savefig('C:/Users/Brendan/Desktop/1600_spectra_resids_slopes_by_freq_log_M1.pdf', format='pdf', bbox_inches='tight')
"""



#resids_full = np.zeros((149))

## load slopes of residual curves by frequency (if loading from file)
#UV_slopes = np.load('C:/Users/Brendan/Desktop/UV_calibration_slopes.npy')

UV_calibrated_params = np.zeros((2,7,bstrp_param_arr.shape[2]))

UV_calibrated_params[0,6] = jul_dates
UV_calibrated_params[1,6] = jul_dates

## loop through wavelengths 
for w in range(2):
    
    ## loop through dates
    for i in range(bstrp_param_arr.shape[2]):
    #for i in range(26):
        
        A, n, C, P, fp, fw = bstrp_param_arr[w,:6,i]  # extract averaged parameters for each dataset
    
        fp = np.log((1./(fp*60.)))  # change Gaussian location to frequency from period
        
        m2_fit = GaussPowerBase(freqs, A,n,C,P,fp,fw)  # generate averaged spectra  
        #m2P_fit = PowerLaw(freqs, A, n, C)
        #m2G_fit = Gauss(freqs, P, fp, fw)
        
        t = [(jul_dates[i] - jul_dates[0]) for q in range(149)]  # determine time between given date and baseline
        diff = UV_slopes[w]*t  # calculate calibration to be applied
        
        f = freqs
        m2_fit_rev = m2_fit - diff  # apply calibration curve to averaged spectra
        
        ## fit model M2 to the calibrated curve
        # assign equal weights to all parts of the curve
        df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
        df2 = np.zeros_like(f)
        df2[0:len(df)] = df
        df2[len(df2)-1] = df2[len(df2)-2]
        ds = df2
        
        M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
        M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
        
        nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, freqs, m2_fit_rev, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp2  # unpack calibrated fitting parameters     
        fit2 = GaussPowerBase(freqs, A2, n2, C2, P2, fp2, fw2)
        
        UV_calibrated_params[w,0:6,i] = nlfit_gp2
        
        #print n, (1./np.exp(fp))/60.
        #print n2, (1./np.exp(fp2))/60.
        
        #resids_temp = m2_fit - m2P
        
        #resids_full = np.vstack((resids_full, resids_temp))
        
        
        """
        ## plot raw and calibrated spectra
        #fig = plt.figure(figsize=(12,10))
        #plt.title(r'1600 $\AA$', fontsize=23, y=1.01)
        plt.loglog(freqs, m2_fit, linewidth=2., label='%s' % greg_dates[i])
        plt.loglog(freqs, m2_fit_rev, linewidth=2.,linestyle='dashed', label='Calibrated')
        #plt.loglog(freqs, fit2, linewidth=2.)
        #plt.loglog(freqs, m2P, 'r--', linewidth=2.)
        plt.ylim((10**-4.7,10**-0.5))
        plt.xlim((10**-4.,10**-1.5))
        #plt.ylim((10**-3.,10**-1.))
        #plt.xlim((10**-3.,10**-2.))
        plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', linewidth=3.)
        plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', linewidth=3.)
        plt.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
        #plt.savefig('C:/Users/Brendan/Desktop/1600_spectra_calibrated.pdf', format='pdf', bbox_inches='tight')
        """

#np.save('C:/Users/Brendan/Desktop/UV_calibrated_params_Nov_M2.npy', UV_calibrated_params)