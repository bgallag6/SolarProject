# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 12:12:45 2017

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
    
spectra_array = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/193/20130626_193_-450_-200i_-200_200j_spectra.npy')
## load in array of segment-averaged pixel FFTs
SPECTRA = spectra_array

print "The region size is %ii x %ij" % (SPECTRA.shape[0], SPECTRA.shape[1])
print "%i frequencies were evaluated in the FFT" % SPECTRA.shape[2] 

num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used
    
# determine frequency values that FFT will evaluate
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # add as argument, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]


start = timer()
T1 = 0


#for l in range(0,SPECTRA.shape[0]):
for l in range(145,146):
    
    #for m in range(0,SPECTRA.shape[1]):
    for m in range(190,191):
        
                                        
        f = freqs  # frequencies
        s = spectra_array[l][m]  # fourier power
        
       
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
            nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
           
        
        except RuntimeError:
            print("Error M1 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M1 - inf/NaN - %i, %i" % (l,m))

      
        A, n, C = nlfit_l  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
        
        
        ## possibly use previous pixel's parameters as initial guesses for current pixel (issues creating wierd banding in images)
        
        ## fit data to combined power law plus gaussian component model
        #"""        
        try:                                 
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
            #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
            
            # change method to 'dogbox' and increase max number of function evaluations to 3000
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.43], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            
        except RuntimeError:
            print("Error M2 - curve_fit failed - %i, %i" % (l,m))
        
        except ValueError:
            print("Error M2 - inf/NaN - %i, %i" % (l,m))
        #"""
        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
        """
        #m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]
        """
        uncertainties = dA2, dn2, dC2, dP2, dfp2, dfw2  # do we want to keep a global array of uncertainties?  
        
        # create model functions from fitted parameters
        #m1_fit = PowerLaw(f_fit, A, n, C)
        m1_fit = PowerLaw(f, A, n, C)
        amp_scale = PowerLaw(np.exp(fp2), A, n, C)  # to extract the gaussian-amplitude scaling factor
        #m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
        m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
        #s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
        m2P_fit = PowerLaw(f, A2, n2, C2)  # only need if plotting
        m2G_fit = Gauss(f, P2, fp2, fw2)  # only need if plotting
        
        #diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
        #diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
        
        #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
        nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
        m2_param = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
        m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22)     
        residsM22 = (s - m2_fit2)
        chisqrM22 = ((residsM22/ds)**2).sum()
        redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6)  
        
        m2P_fit = PowerLaw(f, A22, n22, C22)  # only need if plotting
        m2G_fit = Gauss(f, P22, fp22, fw22)  # only need if plotting           
        
        #residsM2 = (s - s_fit_gp_full)
        residsM2 = (s - m2_fit)
        chisqrM2 = ((residsM2/ds)**2).sum()
        redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
        
        residsM1 = (s - m1_fit)
        chisqrM1 =  ((residsM1/ds)**2).sum()
        redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       
        
        #f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))        
        f_test = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))   
        # Plot models + display combined-model parameters + uncertainties
        #"""
        fig = plt.figure(figsize=(20,15))
        plt.title('Index Bounds [0.3, 6.0], dogbox + trf: Pixel %ii, %ij' % (l,m), y = 1.01, fontsize=25)
        plt.ylim((10**-8,10**1))
        plt.xlim((10**-5,10**-1))
        plt.loglog(f,s,'k')
        plt.loglog(f, m1_fit, label='Power Law - M1')
        plt.loglog(f, m2P_fit, 'g', label='Power Law - M2')
        plt.loglog(f, m2G_fit, 'g--', label='Gaussian - M2')
        #plt.loglog(f, m2_fit, 'r', label='Combined - M2')
        plt.loglog(f, m2_fit2, 'purple', label='Combined - M22')
        plt.xlabel('Frequency (Hz)', fontsize=20, labelpad=10)
        plt.ylabel('Power', fontsize=20, labelpad=10)
        plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
        plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
        plt.text(0.01, 10**0., 'A = {0:0.3e}$\pm${1:0.3e}'.format(m2_param[0], uncertainties[0]), fontsize=15)
        plt.text(0.01, 10**-0.3, 'n = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[1], uncertainties[1]), fontsize=15)
        plt.text(0.01, 10**-0.6, 'C = {0:0.3e}$\pm${1:0.3e}'.format(m2_param[2], uncertainties[2]), fontsize=15)
        plt.text(0.01, 10**-0.9, 'P = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[3], uncertainties[3]), fontsize=15)
        plt.text(0.01, 10**-1.2, 'fp = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[4], uncertainties[4]), fontsize=15)
        plt.text(0.01, 10**-1.5, 'fw = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[5], uncertainties[5]), fontsize=15)
        plt.text(0.01, 10**-1.8, 'f_test = {0:0.3f}'.format(f_test), fontsize=15)
        plt.text(0.01, 10**-2.1, 'chi2 = {0:0.4f}'.format(redchisqrM2), fontsize=15)
        plt.text(0.01, 10**-2.4, 'chi22 = {0:0.4f}'.format(redchisqrM22), fontsize=15)
        plt.legend(loc='upper left', prop={'size':15})
        plt.show()
        #plt.savefig('C:/Users/Brendan/Desktop/bounds_test/%ii_%ij_6_double_optimize.jpeg' % (l,m))
        #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
        #plt.close()
        #"""
        
     