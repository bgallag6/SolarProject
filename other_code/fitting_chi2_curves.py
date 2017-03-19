# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 10:43:07 2017

@author: Brendan
"""

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
        
        
        A2p = 6.65e-09
        n2p = 2.053
        C2p = 2.643e-04
        P2p = 3.146042992254266210e-04
        fp2p = -6.446326188774092358e+00
        fw2p = 0.8        
        
        a_step = A2p/100. 
        n_step = n2p/100. 
        c_step = C2p/100. 
        p_step = P2p/100. 
        fp_step = fp2p/100. 
        fw_step = fw2p/100. 
        
        
        #low_chi = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1]))
        #low_chi = np.zeros((np.random_integers(20,30),np.random_integers(20,30)))
        low_chi = 50
        countT = 0
        

        redchisqrgp = np.zeros((6,100))
        n = np.zeros((6,100))
        
        A2 = A2p
        C2 = C2p
        n2 = n2p
        P2 = P2p
        fp2 = fp2p
        fw2 = fw2p
                
        
        params_step = [a_step, n_step, c_step, p_step, fp_step, fw_step]
        
        for j in range(6):
            
            for i in range(-50,50):
                
                params = [A2, n2, C2, P2, fp2, fw2]                
                params[j] = params[j] + (i*params_step[j])
                               
                #n2 = 2.0 + (i*n_step)
                
                s_fit_gp_full = GaussPowerBase(f, params[0],params[1],params[2],params[3],params[4],params[5])
                residsgp = (s - s_fit_gp_full)
                redchisqrgp[j][i+50] = ((residsgp/ds)**2).sum()/float(f.size-6)                                                             
                
                n[j][i+50] = params[j]
        

            
            fig = plt.figure(figsize=(20,15))
            #plt.loglog(f,s,'k')
            plt.plot(n[j],redchisqrgp[j])
                  
fig = plt.figure(figsize=(20,15))
#plt.loglog(f,s,'k')
plt.plot(n[1],redchisqrgp[1])
plt.ylim(0,0.5)                              
                                
             