# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 18:01:04 2016

@author: Brendan
"""

"""
# notes

# update 1/10 : took out a bunch of unneeded modules
#
# added array for uncertainties and M1M2 model difference
# 
# changed FFT_6segments variable name to spec_arr
#
# changed x_interp / f_fit to = freqs, rather than interpolate to 2x or more points
# also got rid of making above nparray, since freqs is already
# (need a run through to see what heatmaps look like?)
# run through with 1600 rebin4 was identical - except for power law index histogram shifted to the right, but perfectly
# (try with one of jacks full regions)
# tried with 20120923 171 full region - looks great

# changed M1 model fit method to "dogbox" - looks like saves around 15% time

# changed M1 bounds from bounds=([-4.34,0.5,-8.68], [2.,6.,2.])
# to ([-0.1, 0.5, -0.1], [[0.1, 4., 0.1])  ** saved around 15% time - (actually maybe not)

# replaced #'s of bounds within curve_fit call with arrays

# update 1/12: added estimated time remaining

# t_interp might be off by one point in everything, not sure if matters
# wow, t_interp actually not needed at all - thats embarassing
# took out completely
# took out TIME array load completely

# possible arguments to function call = SPECTRA, TIME, pixel-box size?, custom parameters?

# put in module 12:40 AM 1/13
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
from scipy import fftpack
from statsmodels.nonparametric.smoothers_lowess import lowess
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

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)
    

## load in array of segment-averaged pixel FFTs
SPECTRA = np.load('C:/Users/Brendan/Desktop/SDO/spectra_array_FFT6_20130530_1600_2300_2600i_2200_3000j_rebin4_rev_t_interp.npy')
#SPECTRA = np.load('F:/Users/Brendan/Desktop/SolarProject/spectra_array_FFT/spectra_array_FFT6_20130530_1600_2300_2600i_2200_3000j_rebin4_no_median.npy')

## load in time array
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/time_arrays/SDO_20130530_1600A_2300_2600i_2200_3000j_float_time.npy')
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/time_arrays/SDO_20120923_171A_(528)_(132)x_(100)_100y_time.npy')


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


# initialize arrays to hold temporary results for calculating geometric average
p_geometric = np.zeros((num_freq))  # would pre-allocating help? (seems to)
temp = np.zeros((9,num_freq))  # maybe have 3x3 to be generalized

# initialize arrays to hold parameter values, also each pixel's spectra and combined model fit - for tool
diffM1M2 = np.zeros((SPECTRA.shape[0]-2,SPECTRA.shape[1]-2))  # dont really use - get rid of?
params = np.zeros((7,SPECTRA.shape[0]-2,SPECTRA.shape[1]-2))
spectra = np.zeros((SPECTRA.shape[0]-2,SPECTRA.shape[1]-2,SPECTRA.shape[2]))
#M2_fit = np.zeros((SPECTRA.shape[0]-2,SPECTRA.shape[1]-2,(len(freqs)+1)/2))  # would save storage / memory space
M2_fit = np.zeros((SPECTRA.shape[0]-2,SPECTRA.shape[1]-2,SPECTRA.shape[2]))

Uncertainties = np.zeros((6,SPECTRA.shape[0]-2,SPECTRA.shape[1]-2))

#visual_avg = np.average(data,axis=0)  # add visual image to params array - easy access for heatmaps / tool
#params[7] = visual_avg

start = timer()
T1 = 0


### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

for l in range(1,SPECTRA.shape[0]-1):
#for l in range(25,26):
    #print l
    for m in range(1,SPECTRA.shape[1]-1):
    #for m in range(80,85):
        
        temp[0] = np.log10(SPECTRA[l-1][m-1])
        temp[1] = np.log10(SPECTRA[l-1][m])
        temp[2] = np.log10(SPECTRA[l-1][m+1])
        temp[3] = np.log10(SPECTRA[l][m-1])
        temp[4] = np.log10(SPECTRA[l][m])
        temp[5] = np.log10(SPECTRA[l][m+1])
        temp[6] = np.log10(SPECTRA[l+1][m-1])
        temp[7] = np.log10(SPECTRA[l+1][m])
        temp[8] = np.log10(SPECTRA[l+1][m+1])

        temp2 = np.sum(temp, axis=0)
        p_geometric = temp2 / 9.
        p_geometric = np.power(10,p_geometric)
                            
        
        f = freqs  # frequencies
        s = p_geometric  # fourier power
        
        #ds = (1./f**2.2)/1000
        ds = p_geometric*0.1  # set the error / variance estimate to a constant percentage of the spectra power-values
        
        # create points to fit model with final parameters 
        #f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],(len(freqs)+1)/2)  # would save storage / memory space
        f_fit = freqs       
    
                
                
        ### fit data to models using SciPy's Levenberg-Marquart method
        
        try:
            # initial guesses for fitting parameters
            #P0 = [0.000, 2.0, 0.00003, 0.0022, -6.5, 0.5]
            M1_low = [-0.1, 0.5, -0.1]
            M1_high = [0.1, 4., 0.1]
            nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
           
        
        except RuntimeError:
            print("Error M1 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M1 - inf/NaN - %i, %i" % (l,m))

      
        A, n, C = nlfit_l  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
        
        
        ## possibly use previous pixel's parameters as initial guesses for current pixel (issues creating wierd banding in images)
                        
        """
        if m > 1:
            P0 = [pl_A[l-1][m-2], slopes[l-1][m-2], pl_C[l-1][m-2], gauss[l-1][m-2], gauss_loc[l-1][m-2], gauss_wid[l-1][m-2]]
            #P0 = [0.0, 1.02, 0.001, 0.001, -4.68, 0.79]
            try:
                #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
        
            except RuntimeError:
                print("Error M2 - curve_fit failed")
            
            except ValueError:
                print("Error M2 - inf/NaN - %i, %i" % (l,m))
                
        elif m == 1 and l > 1:
            P0 = [pl_A[l-2][m-1], slopes[l-2][m-1], pl_C[l-2][m-1], gauss[l-2][m-1], gauss_loc[l-2][m-1], gauss_wid[l-2][m-1]]
            try:
                #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
        
            except RuntimeError:
                print("Error M2 - curve_fit failed")
            
            except ValueError:
                print("Error M2 - inf/NaN - %i, %i" % (l,m))
                
        else:        
        #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0=[nlfit_l[0], nlfit_l[1], nlfit_l[2], nlfit_g[0], nlfit_g[1], nlfit_g[2]], sigma=ds)
            try:
                #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
        
            except RuntimeError:
                print("Error M2 - curve_fit failed")
            
            except ValueError:
                print("Error M2 - inf/NaN - %i, %i" % (l,m))
        """     
        
        ## fit data to combined power law plus gaussian component model
        #"""        
        try:
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)                                  
            M2_low = [-0.1, 0.5, -0.1, 0.00001, -6.5, 0.05]
            M2_high = [0.1, 4., 0.1, 0.2, -4.6, 0.8]
            # change method to 'dogbox' and increase max number of function evaluations to 3000
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            
        except RuntimeError:
            print("Error M2 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M2 - inf/NaN - %i, %i" % (l,m))
        #"""
        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
        m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]
        uncertainties = dA2, dn2, dC2, dP2, dfp2, dfw2  # do we want to keep a global array of uncertainties?  
        
        
        uncertainties_arr = [dA2, dn2, dC2, dP2, dfp2, dfw2]  # not sure if want to keep these
        Uncertainties[:, l-1, m-1] = uncertainties_arr
        
        
        # create model functions from fitted parameters
        m1_fit = PowerLaw(f_fit, A, n, C)
        m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
        s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
        m2P_fit = PowerLaw(f_fit, A2, n2, C2)
        m2G_fit = Gauss(f_fit, P2, fp2, fw2)
        
        diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
        diffM1M2[l-1][m-1] = np.sum(diffM1M2_temp)  # sum of squared differences 
        
        residsgp = (s - s_fit_gp_full)
        redchisqrgp = ((residsgp/ds)**2).sum()/float(f.size-6)
        
        # populate array with parameters
        params[0][l-1][m-1] = A2
        params[1][l-1][m-1] = n2
        params[2][l-1][m-1] = C2
        params[3][l-1][m-1] = P2
        params[4][l-1][m-1] = fp2
        params[5][l-1][m-1] = fw2
        params[6][l-1][m-1] = redchisqrgp
        
        # populate arrays holding spectra / fits
        spectra[l-1][m-1] = p_geometric
        M2_fit[l-1][m-1] = m2_fit
        
        
        
        # Plot models + display combined-model parameters + uncertainties
        """
        fig = plt.figure(figsize=(20,15))
        plt.title('SDO AIA 304.0 Angstrom 20120923 - 3 Segments, 3x3 Pixel-Box Averaging 598 interp', y = 1.01, fontsize=25)
        plt.ylim((10**-8,10**1))
        plt.xlim((10**-5,10**-1))
        plt.loglog(f,s,'k')
        plt.loglog(f_fit, m1_fit, label='Power Law - M1')
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
    
    # estimate time remaining and print to screen   
    T = timer() - start
    T2 = T - start - T1
    if l == 1:
        T_est = T2*(SPECTRA.shape[0]-2)    
    T_est2 = T2*((SPECTRA.shape[0]-2)-l+1)
    print "Currently on row %i of %i, estimated time remaining: %i seconds" % (l, SPECTRA.shape[0]-2, T_est2)
    T1 = T

# print estimated and total program time to screen        
print "Beginning Estimated time = %i sec" % T_est
T_act = timer() - start
print "Actual total time = %i sec" % T_act           

