# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 20:51:12 2017

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
import h5py
from scipy import fftpack  # doesnt work in module when called here???
import matplotlib.pylab as plt
from astropy.convolution import convolve, Box1DKernel
from matplotlib import cm
from numpy.random import randn
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.ticker import LogFormatterMathtext
from timeit import default_timer as timer
from scipy.stats import f
import matplotlib.patches as patches
from sunpy.map import Map
from scipy.stats import f as ff

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
    
#spectra_array = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/193/20130626_193_-450_-200i_-200_200j_spectra.npy')
#spectra_array = np.load('C:/Users/Brendan/Desktop/project_files/20130626_171_-500_500i_-500_600j_spectra_arth.npy')
#visual = np.load('C:/Users/Brendan/Desktop/solar_final/20130626_171_-500_500i_-500_600j_visual.npy')
spec_171 = np.load('C:/Users/Brendan/Desktop/spec_slice/171_323.npy')
spec_304 = np.load('C:/Users/Brendan/Desktop/spec_slice/211_323.npy')
#spectra_array = np.load('C:/Users/Brendan/Desktop/1600/spectra.npy')
#spectra_array = np.memmap('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20130626/1600_older_spectra_mmap.npy', dtype='float64', mode='r', shape=(1636,1621,299))
#spectra_array = np.memmap('/mnt/data-solar/Gallagher/data_older/20130626/1600_other/spectra_mmap.npy', dtype='float64', mode='r', shape=(1636,1621,299))
#visual = np.load('C:/Users/Brendan/Desktop/1600/visual_1600.npy')
#param = np.load('C:/Users/Brendan/Desktop/1600/param_1600.npy')
param171 = np.load('C:/Users/Brendan/Desktop/20130626_final/20130626_171_-500_500i_-500_600j_param_slope6_arthm.npy')
param304 = np.load('C:/Users/Brendan/Desktop/20130626_final/20130626_211_-500_500i_-500_600j_param_slope6_arthm.npy')
p_loc_171 = param171[4]
p_loc_304 = param304[4]
#p_loc = 1./np.exp(param[4])

#vis = visual[0]

# generate p-value heatmap
df1, df2 = 3, 6
p_val171 = ff.sf(param171[6], df1, df2)
p_val304 = ff.sf(param304[6], df1, df2)


mask_thresh = 0.005    
   
loc_mask171 = np.copy(p_loc_171)
loc_mask304 = np.copy(p_loc_304)

#h_min_amp = np.percentile(h_map[3],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
#h_max_amp = np.percentile(h_map[3],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)

#count = 0

for i in range(p_val171.shape[0]):
        for j in range(p_val171.shape[1]):
            if p_val171[i][j] > mask_thresh:
                #count += 1
                loc_mask171[i][j] = np.NaN
p_loc_171 = loc_mask171[125:525,525:950]
p_loc_171 = 1./np.exp(p_loc_171)

for i in range(p_val304.shape[0]):
        for j in range(p_val304.shape[1]):
            if p_val304[i][j] > mask_thresh:
                #count += 1
                loc_mask304[i][j] = np.NaN
p_loc_304 = loc_mask304[125:525,525:950]
p_loc_304 = 1./np.exp(p_loc_304)

## load in array of segment-averaged pixel FFTs

num_freq = 299  # determine nubmer of frequencies that are used
    
# determine frequency values that FFT will evaluate
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # add as argument, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]


start = timer()
T1 = 0

#m2 = [790, 767, 757, 765, 867, 525, 762, 649, 722, 485, 743, 708, 744, 14] #use
#l2 = [235, 547, 319, 325, 864, 551, 1529, 1548, 1441, 535, 322, 352, 272, 330] #use

#m2 = [722]
#l2 = [1441]

#m2 = [14, 58, 283, 512, 629, 810, 909, 901, 1342, 767, 873]
#l2 = [330, 394, 482, 538, 379, 293, 400, 409, 1239, 1160, 1261]

#m2 = [743, 708, 525, 757, 765, 722, 867]
#l2 = [322, 352, 551, 319, 325, 1441, 864]


y_rang = [321,322]
#x_rang = [685,800]
x_rang = [723,724]

#spectra_points = np.zeros((3,299))

#for l in range(1278,1279):
#for l in range(y_rang[0],y_rang[1]):
for l in range(1):
#for l in range(1):
    
    for m in range(x_rang[0],x_rang[1]):
    #for m in range(685,686):
    #for m in range(900,901):
    #for m in range(len(m2)):    
        
                                        
        f = freqs  # frequencies
        s171 = spec_171[m-x_rang[0]+185]  # fourier power
        #s = spectra_array[l2[m]][m2[m]]  # fourier power
        #spectra_points[m-x_rang[0]] = s
       
        # assign equal weights to all parts of the curve
        df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
        df2 = np.zeros_like(f)
        df2[0:len(df)] = df
        df2[len(df2)-1] = df2[len(df2)-2]
        ds = df2
        #ds = 0.1*s
        
                                               
        ### fit data to models using SciPy's Levenberg-Marquart method
        
        try:
            # initial guesses for fitting parameters
            M1_low = [-0.002, 0.3, -0.01]
            M1_high = [0.002, 6., 0.01]
            nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s171, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
           
        
        except RuntimeError:
            print("Error M1 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M1 - inf/NaN - %i, %i" % (l,m))

      
        A, n, C = nlfit_l  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
        
        
        ## possibly use previous pixel's parameters as initial guesses for current pixel (issues creating wierd banding in images)
        
        ## fit data to combined power law plus gaussian component model
               
        try:                                 
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
            #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
            
            # change method to 'dogbox' and increase max number of function evaluations to 3000
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.43], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s171, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            
        except RuntimeError:
            print("Error M2 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M2 - inf/NaN - %i, %i" % (l,m))
        
        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
        
          
        
        # create model functions from fitted parameters
        #m1_fit = PowerLaw(f_fit, A, n, C)
        m1_fit171 = PowerLaw(f, A, n, C)
        amp_scale = PowerLaw(np.exp(fp2), A, n, C)  # to extract the gaussian-amplitude scaling factor
        #m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
        m2_fit171 = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
        #s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
        #m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]

        
        #diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
        #diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
        
        
        nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s171, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
        #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
        dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        m2_param171 = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
        m2_fit2171 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22) 
        uncertainties = dA22, dn22, dC22, dP22, dfp22, dfw22  # do we want to keep a global array of uncertainties?
        
        residsM22 = (s171 - m2_fit2171)
        chisqrM22 = ((residsM22/ds)**2).sum()
        redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6)  
        
        m2P_fit171 = PowerLaw(f, A22, n22, C22)  # only need if plotting
        m2G_fit171 = Gauss(f, P22, fp22, fw22)  # only need if plotting           
              
        
        #residsM2 = (s - s_fit_gp_full)
        residsM2 = (s171 - m2_fit171)
        chisqrM2 = ((residsM2/ds)**2).sum()
        redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
        
        residsM1 = (s171 - m1_fit171)
        chisqrM1 =  ((residsM1/ds)**2).sum()
        redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       
        
        #f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))        
        f_test171 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))   
        # Plot models + display combined-model parameters + uncertainties
        
        
        s304 = spec_304[m-x_rang[0]+185]  # fourier power
        #s = spectra_array[l2[m]][m2[m]]  # fourier power
        #spectra_points[m-x_rang[0]] = s
       
                                              
        ### fit data to models using SciPy's Levenberg-Marquart method
        
        try:
            # initial guesses for fitting parameters
            M1_low = [-0.002, 0.3, -0.01]
            M1_high = [0.002, 6., 0.01]
            nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s304, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
           
        
        except RuntimeError:
            print("Error M1 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M1 - inf/NaN - %i, %i" % (l,m))

      
        A, n, C = nlfit_l  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
        
        
        ## possibly use previous pixel's parameters as initial guesses for current pixel (issues creating wierd banding in images)
        
        ## fit data to combined power law plus gaussian component model
               
        try:                                 
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
            #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
            
            # change method to 'dogbox' and increase max number of function evaluations to 3000
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.43], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s304, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            
        except RuntimeError:
            print("Error M2 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M2 - inf/NaN - %i, %i" % (l,m))
        
        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
        
          
        
        # create model functions from fitted parameters
        #m1_fit = PowerLaw(f_fit, A, n, C)
        m1_fit304 = PowerLaw(f, A, n, C)
        amp_scale = PowerLaw(np.exp(fp2), A, n, C)  # to extract the gaussian-amplitude scaling factor
        #m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
        m2_fit304 = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
        #s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
        #m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]

        
        #diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
        #diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
        
        
        nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s304, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
        #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
        dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        m2_param304 = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
        m2_fit2304 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22) 
        uncertainties = dA22, dn22, dC22, dP22, dfp22, dfw22  # do we want to keep a global array of uncertainties?
        
        residsM22 = (s304 - m2_fit2304)
        chisqrM22 = ((residsM22/ds)**2).sum()
        redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6)  
        
        m2P_fit304 = PowerLaw(f, A22, n22, C22)  # only need if plotting
        m2G_fit304 = Gauss(f, P22, fp22, fw22)  # only need if plotting           
              
        
        #residsM2 = (s - s_fit_gp_full)
        residsM2 = (s304 - m2_fit304)
        chisqrM2 = ((residsM2/ds)**2).sum()
        redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
        
        residsM1 = (s304 - m1_fit304)
        chisqrM1 =  ((residsM1/ds)**2).sum()
        redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       
        
        #f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))        
        f_test304 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))   
        # Plot models + display combined-model parameters + uncertainties
        
        fig = plt.figure(figsize=(16,25))
        ax1 = plt.subplot2grid((2,16),(0, 0), colspan=8, rowspan=1)
        #plt.title('Power-Law Dominated : Pixel %ii, %ij' % (l2[m],m2[m]), y = 1.01, fontsize=25)
        ax1.set_title('171A: Pixel %ix, %iy' % (m,321), y = 1.01, fontsize=21)
        #ax1.set_title('1600A: Pixel %ix, %iy' % (m,l), y = 1.01, fontsize=21)
        ax1.set_ylim((10**-5,10**0))
        ax1.set_xlim((10**-4.1,10**-1.2))
        ax1.loglog(f,s171,'k')
        ax1.loglog(f, m1_fit171, label='M1 - Power Law')
        ax1.loglog(f, m2P_fit171, 'g', label='M2 - Power Law')
        ax1.loglog(f, m2G_fit171, 'g--', label='M2 - Gaussian')
        #plt.loglog(f, m2_fit, 'r', label='Combined - M2')
        ax1.loglog(f, m2_fit2171, 'purple', label='M2 - Combined')
        #ax1.set_xlabel('Frequency (Hz)', fontsize=13, labelpad=7)
        ax1.set_ylabel('Power', fontsize=13, labelpad=7)
        ax1.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
        ax1.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
        
        ax1.text(0.008, 10**-0.45, r'$A$ = {0:0.2e}'.format(m2_param171[0]), fontsize=15)
        ax1.text(0.008, 10**-0.70, r'$n$ = {0:0.3f}'.format(m2_param171[1]), fontsize=15)
        ax1.text(0.008, 10**-0.95, r'$R$ = {0:1.0f} '.format((m2_param171[2]/m2_param171[0])**(-1./m2_param171[1])), fontsize=15)
        ax1.text(0.008, 10**-1.2, r'$\alpha$ = {0:0.2e}'.format(m2_param171[3]), fontsize=15)
        ax1.text(0.008, 10**-1.45, r'$\beta$ = {0:1.0f} [s]'.format(1./np.exp(m2_param171[4])), fontsize=15)
        ax1.text(0.008, 10**-1.70, r'$\sigma$ = {0:0.3f}'.format(m2_param171[5]), fontsize=15)
        #plt.text(0.01, 10**-2.4, r'$\chi^2$: Dogbox + trf = {0:0.4f}'.format(redchisqrM22), fontsize=15)
        #ax1.legend(loc='upper left', prop={'size':15})
        ax1.legend(loc='lower left', prop={'size':12})

        #""" Gaussian Location Side-Plot
        ax2 = plt.subplot2grid((2,16),(0, 8), colspan=8, rowspan=1) 
        ax2.imshow(np.flipud(p_loc_171), cmap='jet_r')
        #ax2.set_title('171A: Visual Average', y = 1.01, fontsize=17)
        ax2.set_title(r'171$\AA$: Gaussian Location - Masked', y = 1.01, fontsize=21)
        rect = patches.Rectangle((x_rang[0]-525,y_rang[0]-125), (x_rang[1]-x_rang[0]), 7, color='white', fill=False, linewidth=2)
        rect3 = patches.Rectangle(((m-3-525),y_rang[0]-125), 7, 7, color='red', fill=True)
        rect2 = patches.Rectangle(((m-2-525),y_rang[0]-124), 5, 5, color='white', fill=True)  
        ax2.add_patch(rect)
        ax2.add_patch(rect3)
        ax2.add_patch(rect2)
        #"""
        
        ax3 = plt.subplot2grid((2,16),(1, 0), colspan=8, rowspan=1)
        #plt.title('Power-Law Dominated : Pixel %ii, %ij' % (l2[m],m2[m]), y = 1.01, fontsize=25)
        ax3.set_title('211A: Pixel %ix, %iy' % (m,321), y = 1.01, fontsize=21)
        #ax1.set_title('1600A: Pixel %ix, %iy' % (m,l), y = 1.01, fontsize=21)
        ax3.set_ylim((10**-5,10**0))
        ax3.set_xlim((10**-4.1,10**-1.2))
        ax3.loglog(f,s304,'k')
        ax3.loglog(f, m1_fit304, label='M1 - Power Law')
        ax3.loglog(f, m2P_fit304, 'g', label='M2 - Power Law')
        ax3.loglog(f, m2G_fit304, 'g--', label='M2 - Gaussian')
        #plt.loglog(f, m2_fit, 'r', label='Combined - M2')
        ax3.loglog(f, m2_fit2304, 'purple', label='M2 - Combined')
        ax3.set_xlabel('Frequency (Hz)', fontsize=13, labelpad=7)
        ax3.set_ylabel('Power', fontsize=13, labelpad=7)
        ax3.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
        ax3.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
        
        ax3.text(0.008, 10**-0.45, r'$A$ = {0:0.2e}'.format(m2_param304[0]), fontsize=15)
        ax3.text(0.008, 10**-0.70, r'$n$ = {0:0.3f}'.format(m2_param304[1]), fontsize=15)
        ax3.text(0.008, 10**-0.95, r'$R$ = {0:1.0f} '.format((m2_param304[2]/m2_param304[0])**(-1./m2_param304[1])), fontsize=15)
        ax3.text(0.008, 10**-1.20, r'$\alpha$ = {0:0.2e}'.format(m2_param304[3]), fontsize=15)
        ax3.text(0.008, 10**-1.45, r'$\beta$ = {0:1.0f} [s]'.format(1./np.exp(m2_param304[4])), fontsize=15)
        ax3.text(0.008, 10**-1.70, r'$\sigma$ = {0:0.3f}'.format(m2_param304[5]), fontsize=15)
        #plt.text(0.01, 10**-2.4, r'$\chi^2$: Dogbox + trf = {0:0.4f}'.format(redchisqrM22), fontsize=15)
        #ax1.legend(loc='upper left', prop={'size':15})
        ax3.legend(loc='lower left', prop={'size':12})
        
        """ Visual Image Side-Plot
        ax2 = plt.subplot2grid((1,26),(0, 14), colspan=12, rowspan=1)
        v_min = np.percentile(visual[0],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        v_max = np.percentile(visual[0],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)  
        #ax2.imshow(visual[0], cmap='sdoaia171', vmin=v_min, vmax=v_max)
        ax2.imshow(visual[0], cmap='sdoaia1600', vmin=v_min, vmax=v_max)
        #ax2.set_title('171A: Visual Average', y = 1.01, fontsize=17)
        ax2.set_title('1600A: Visual Average', y = 1.01, fontsize=21)
        ax2.set_xlim(0, visual[0].shape[1])
        ax2.set_ylim(0, visual[0].shape[0])
        rect = patches.Rectangle((x_rang[0],y_rang[0]), (x_rang[1]-x_rang[0]), 9, color='white', fill=False, linewidth=2)
        rect2 = patches.Rectangle(((m-3),y_rang[0]), 9, 9, color='red', fill=True)
        ax2.add_patch(rect)
        ax2.add_patch(rect2)
        """
        
        #""" Gaussian Location Side-Plot
        ax4 = plt.subplot2grid((2,16),(1, 8), colspan=8, rowspan=1) 
        ax4.imshow(np.flipud(p_loc_304), cmap='jet_r')
        #ax2.set_title('171A: Visual Average', y = 1.01, fontsize=17)
        ax4.set_title(r'211$\AA$: Gaussian Location - Masked', y = 1.01, fontsize=21)
        rect = patches.Rectangle((x_rang[0]-525,y_rang[0]-125), (x_rang[1]-x_rang[0]), 7, color='white', fill=False, linewidth=2)
        rect3 = patches.Rectangle(((m-3-525),y_rang[0]-125), 7, 7, color='red', fill=True)
        rect2 = patches.Rectangle(((m-2-525),y_rang[0]-124), 5, 5, color='white', fill=True)  
        ax4.add_patch(rect)
        ax4.add_patch(rect3)
        ax4.add_patch(rect2)
        
        #plt.savefig('C:/Users/Brendan/Desktop/spectra_points/193_%ii_%ij.pdf' % (l2[m],m2[m]), format='pdf')
        #np.save('/mnt/data-solar/Gallagher/DATA/1600_test/1600_prev_spectra.npy', spectra_points)
        #plt.savefig('C:/Users/Brendan/Desktop/171_211_slice/171_%ix_%iy.jpeg' % (m,l))
        #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
        #plt.close()
#"""