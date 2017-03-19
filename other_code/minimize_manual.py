# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 10:12:15 2017

@author: Brendan
"""

# maybe add, if next iteration chi is larger than previous, break

# maybe do a comprehensive fit on initial pixel - then use seed parameters
# and search in range around it

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from matplotlib.widgets import  RectangleSelector
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
from scipy import fftpack  # doesnt work in module when called here???
#from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib.pylab as plt
from astropy.convolution import convolve, Box1DKernel
from matplotlib import cm
from numpy.random import randn
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from timeit import default_timer as timer




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
    
SPECTRA = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/193/20130626_193_-450_-200i_-200_200j_spectra.npy')
## load in array of segment-averaged pixel FFTs
spectra_array = SPECTRA

print "The region size is %ii x %ij" % (SPECTRA.shape[0], SPECTRA.shape[1])
print "%i frequencies were evaluated in the FFT" % SPECTRA.shape[2] 

num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used
    
# determine frequency values that FFT will evaluate
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # add as argument, or leave in as constant?
#time_step = 24  # add as argument, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]


# initialize arrays to hold parameter values, also each pixel's combined model fit - for tool
diffM1M2 = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1]))  # dont really use - get rid of?
params = np.zeros((7, SPECTRA.shape[0], SPECTRA.shape[1]))
#M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], (len(freqs)+1)/2))  # would save storage / memory space
M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], SPECTRA.shape[2]))


start = timer()
T1 = 0


### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

#for l in range(0,SPECTRA.shape[0]):
for l in range(145,146):
    #print l
    #for m in range(0,SPECTRA.shape[1]):
    for m in range(190,191):
        
                                        
        f = freqs  # frequencies
        s = spectra_array[l][m]  # fourier power
        
        #ds = (1./f**2.2)/1000
        df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
        df2 = np.zeros_like(f)
        df2[0:len(df)] = df
        df2[len(df2)-1] = df2[len(df2)-2]
        ds = df2        
        #ds = s*0.1  # set the error / variance estimate to a constant percentage of the spectra power-values
        
        # create points to fit model with final parameters 
        #f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],(len(freqs)+1)/2)  # would save storage / memory space
        f_fit = freqs       
        
        M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
        #M2_low = [-0.1, 0.1, -0.1, 0.00001, -6.5, 0.05]  # test on 193 - coronal hole
        M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
        
        param_bin = 2
        
        p1_range = (M2_high[0]-M2_low[0])/10. 
        p2_range = (M2_high[1]-M2_low[1])/10.      
        p3_range = (M2_high[2]-M2_low[2])/10.  
        #p4_range = (M2_high[3]-M2_low[3])/10.      
        #p5_range = (M2_high[4]-M2_low[4])/10.      
        #p6_range = (M2_high[5]-M2_low[5])/10.  
        p4_range = (M2_high[3]-M2_low[3])/100.      
        p5_range = (M2_high[4]-M2_low[4])/100.      
        p6_range = (M2_high[5]-M2_low[5])/100.   
        
        
        n_step = 0.01        
        
        
        #low_chi = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1]))
        #low_chi = np.zeros((np.random_integers(20,30),np.random_integers(20,30)))
        low_chi = 50
        countT = 0
        
        for p1 in range(-5,5):
            for p2 in range(-5,5):
                for p3 in range(-5,5):
                    for p4 in range(-5,5):
                        for p5 in range(-5,5):
                            for p6 in range(-5,5):
                                countT += 1
                                #A2 = M2_low[0] + p1_range*p1
                                #n2 = M2_low[1] + p2_range*p2
                                #C2 = M2_low[2] + p3_range*p3
                                #P2 = M2_low[3] + p4_range*p4
                                #fp2 = M2_low[4] + p5_range*p5
                                #fw2 = M2_low[5] + p6_range*p6
                                
                                A2p = 6.65e-09
                                n2p = 2.053
                                C2p = 2.643e-04
                                P2p = 2.458691276650307417e-04
                                fp2p = -6.191295483972629299e+00
                                fw2p = 0.8
                                
                                A2 = A2p + A2p/10*p1
                                n2 = n2p + n2p/10*p2
                                C2 = C2p + C2p/10*p3
                                P2 = P2p + P2p/10*p4
                                fp2 = fp2p + fp2p/10*p5
                                fw2 = fw2p + fw2p/10*p6
                                                               
                                s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
                                residsgp = (s - s_fit_gp_full)
                                redchisqrgp = ((residsgp/ds)**2).sum()/float(f.size-6)
                                
                                #fig = plt.figure(figsize=(20,15))
                                #plt.loglog(f,s,'k')
                                #plt.loglog(f, s_fit_gp_full, 'b')
                                
                                
                                if redchisqrgp < low_chi:
                                    low_chi = redchisqrgp
                                    A3 = A2
                                    n3 = n2
                                    C3 = C2
                                    P3 = P2
                                    fp3 = fp2
                                    fw3 = fw2
                                    
        
        s_fit_gp_full = GaussPowerBase(f, A3,n3,C3,P3,fp3,fw3)
        fig = plt.figure(figsize=(20,15))
        plt.loglog(f,s,'k')
        plt.loglog(f, s_fit_gp_full, 'b') 
        plt.ylim((10**-8,10**1))
        plt.xlim((10**-5,10**-1))        
        
        
        nlfit_gpT, nlpcov_gpT = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0=[A3,n2,C3,P3,fp3,fw3], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
        A3T, n3T, C3T, P3T, fp3T, fw3T = nlfit_gpT  # unpack fitting parameters
        s_fit_gp_fullT = GaussPowerBase(f, A3T, n3T, C3T, P3T, fp3T, fw3T)
        residsgpT = (s - s_fit_gp_fullT)
        redchisqrgpT = ((residsgpT/ds)**2).sum()/float(f.size-6)
        fig = plt.figure(figsize=(20,15))
        plt.loglog(f,s,'k')
        plt.loglog(f, s_fit_gp_fullT, 'b')   
        plt.ylim((10**-8,10**1))
        plt.xlim((10**-5,10**-1))                     
                                
                                
                
                
        ### fit data to models using SciPy's Levenberg-Marquart method
        ## fit data to combined power law plus gaussian component model
        #"""        
        try:
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)                                  
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            #M2_low = [-0.1, 0.1, -0.1, 0.00001, -6.5, 0.05]  # test on 193 - coronal hole
            M2_high = [0.002, 4., 0.01, 0.2, -4.6, 0.8]
            #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
            # change method to 'dogbox' and increase max number of function evaluations to 3000
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            
        except RuntimeError:
            print("Error M2 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M2 - inf/NaN - %i, %i" % (l,m))
        #"""
        
        A5, n5, C5, P5, fp5, fw5 = nlfit_gp  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
        m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]
        uncertainties = dA2, dn2, dC2, dP2, dfp2, dfw2  # do we want to keep a global array of uncertainties?  
        
        
        uncertainties_arr = [dA2, dn2, dC2, dP2, dfp2, dfw2]  # not sure if want to keep these
        
        
        # create model functions from fitted parameters
        #m1_fit = PowerLaw(f_fit, A, n, C)
        m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
        s_fit_gp_fullT = GaussPowerBase(f, A5,n5,C5,P5,fp5,fw5)
        m2P_fit = PowerLaw(f_fit, A2, n2, C2)
        m2G_fit = Gauss(f_fit, P2, fp2, fw2)
        
        #diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
        #diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
        
        residsgpT = (s - s_fit_gp_fullT)
        redchisqrgpT = ((residsgpT/ds)**2).sum()/float(f.size-6)
        
        # populate array with parameters
        params[0][l][m] = A2
        params[1][l][m] = n2
        params[2][l][m] = C2
        params[3][l][m] = P2
        params[4][l][m] = fp2
        params[5][l][m] = fw2
        #params[6][l][m] = redchisqrgp
        
        # populate array holding model fits
        M2_fit[l][m] = m2_fit
        
        
        
        # Plot models + display combined-model parameters + uncertainties
        """
        fig = plt.figure(figsize=(20,15))
        plt.title('SDO AIA 304.0 Angstrom 20120923 - 3 Segments, 3x3 Pixel-Box Averaging 598 interp', y = 1.01, fontsize=25)
        plt.ylim((10**-8,10**1))
        plt.xlim((10**-5,10**-1))
        plt.loglog(f,s,'k')
        #plt.loglog(f_fit, m1_fit, label='Power Law - M1')
        plt.loglog(f_fit, m2P_fit, 'g', label='Power Law - M2')
        plt.loglog(f_fit, m2G_fit, 'g--', label='Gaussian - M2')
        plt.loglog(f_fit, m2_fit, 'r', label='Combined - M2')
        plt.xlabel('Frequency (Hz)', fontsize=20, labelpad=10)
        plt.ylabel('Power', fontsize=20, labelpad=10)
        plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
        plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
        plt.vlines((1.0/24.),10**-8,10**1, linestyles='solid', label='24 seconds')
        plt.text(0.01, 10**0., 'A = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[0], uncertainties[0]), fontsize=15)
        plt.text(0.01, 10**-0.3, 'n = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[1], uncertainties[1]), fontsize=15)
        plt.text(0.01, 10**-0.6, 'C = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[2], uncertainties[2]), fontsize=15)
        plt.text(0.01, 10**-0.9, 'P = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[3], uncertainties[3]), fontsize=15)
        plt.text(0.01, 10**-1.2, 'fp = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[4], uncertainties[4]), fontsize=15)
        plt.text(0.01, 10**-1.5, 'fw = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[5], uncertainties[5]), fontsize=15)
        plt.legend(loc='upper left', prop={'size':15})
        plt.show()
        #plt.savefig('C:/Users/Brendan/Desktop/PHYS 326/dogbox_test1/20130530_193A_3x3_6seg_%ii_%ij.jpeg' % (l,m))
        #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
        #plt.close()
        """