# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 15:35:17 2017

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
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

spec = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Temp/20120909/1600/spectra.npy')
spec_std = np.load('C:/Users/Brendan/Desktop/spec_std.npy')

# determine frequency values that FFT will evaluate
num_freq = spec.shape[2]  # determine nubmer of frequencies that are used
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # 12 second cadence for the others
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

i = 6

spec3x3 = spec[i:i+3,i:i+3,:]
spec3x3 = np.reshape(spec3x3, (9,1799))
spec3x3_avg = np.average(spec3x3, axis=0)

#spec100x100 = np.reshape(spec, (10000,1799))
#spec100x100_avg = np.log10(spec100x100)
#spec100x100_avg = np.sum(spec100x100_avg, axis=0) / 10000.
#spec100x100_avg = np.power(10,spec100x100_avg)
#spec100x100_avg = np.average(spec100x100, axis=0)

#plt.loglog(freqs, spec[1][1])
#plt.loglog(freqs, spec3x3_avg)
#plt.loglog(freqs, spec100x100_avg)
w = 82
#"""                            
f = freqs  # frequencies
#s = spectra_array[l][m]  # fourier power
s = spec[w][w]  # fourier power
#s = spec3x3_avg  # fourier power
#s = spec100x100_avg  # fourier power

   
# assign equal weights to all parts of the curve
df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
df2 = np.zeros_like(f)
df2[0:len(df)] = df
df2[len(df2)-1] = df2[len(df2)-2]
#ds = df2
#ds = 0.1*s
ds = spec_std

                                       
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
indx = 595      
try:                                 
    M2_low = [-0.002, 0.3, -0.01, 0.000001, -6.5, 0.05]
    M2_high = [0.002, indx/100, 0.01, 0.2, -4.6, 0.8]
    #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
    
    # change method to 'dogbox' and increase max number of function evaluations to 3000
    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.43], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
    
  
except RuntimeError:
    print("Error M1 - curve_fit failed - %i, %i" % (l,m))
        
except ValueError:
    print("Error M1 - inf/NaN - %i, %i" % (l,m))
A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters

# unpack uncertainties in fitting parameters from diagonal of covariance matrix
dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]


  

# create model functions from fitted parameters
m1_fit = PowerLaw(f, A, n, C)
m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)


m2P_fit = PowerLaw(f, A2, n2, C2)  # only need if plotting
m2G_fit = Gauss(f, P2, fp2, fw2)  # only need if plotting
m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]      
uncertainties = dA2, dn2, dC2, dP2, dfp2, dfw2  # do we want to keep a global array of uncertainties?


#diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
#diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 

#M2_low = [0., 0.2, 0., 0., -6.5, 0.05]
#M2_high = [0.002, 6., 0.001, 0.2, -4.6, 0.8]
nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
#nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
m2_param2 = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22) 
uncertainties2 = dA22, dn22, dC22, dP22, dfp22, dfw22  # do we want to keep a global array of uncertainties?

residsM22 = (s - m2_fit2)
chisqrM22 = ((residsM22/ds)**2).sum()
redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6)  

m2P_fit2 = PowerLaw(f, A22, n22, C22)  # only need if plotting
m2G_fit2 = Gauss(f, P22, fp22, fw22)  # only need if plotting           
      

#residsM2 = (s - s_fit_gp_full)
residsM2 = (s - m2_fit)
chisqrM2 = ((residsM2/ds)**2).sum()
redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)

residsM1 = (s - m1_fit)
chisqrM1 =  ((residsM1/ds)**2).sum()
redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       

# Plot models + display combined-model parameters + uncertainties
residsM22 = (s - m2_fit2)
chisqrM22 = ((residsM22/ds)**2).sum()
redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6) 
    
f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))
f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))

  
r_temp = pearsonr(m2_fit2, s)  # calculate r-value correlation coefficient
r = r_temp[0]
#actual = GaussPowerBase(f, 1.3e-12, 2.45, 3.45e-4, 1e-7, -5.5, 0.05)  # only need if plotting

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27
        
fig = plt.figure(figsize=(15,15))
ax = plt.gca()  # get current axis -- to set colorbar 
#plt.title('Power-Law Dominated : Pixel %ii, %ij' % (l2[m],m2[m]), y = 1.01, fontsize=25)
plt.title('Jacks Dataset: Index = 3.0 | No Segment-Averaging | No Pixels Averaged \n Index Upper Bound = %0.2f' % (indx/100.) , y = 1.01, fontsize=25)
plt.ylim((10**-4.7,10**0))
plt.xlim((10**-5.,10**-1.3))
plt.xticks(fontsize=19)
plt.yticks(fontsize=19)
plt.loglog(f,s,'k')
#plt.loglog(f,actual,'red', linewidth=1.3, label='Using 2.45')
#plt.loglog(f, m1_fit, label='M1 - Power Law', linewidth=1.3)
plt.loglog(f, m2P_fit, 'g', label='M2 - Power Law', linewidth=1.3)
plt.loglog(f, m2G_fit, 'g--', label='M2 - Gaussian', linewidth=1.3)
plt.loglog(f, m2_fit, 'purple', label='M2 - Combined', linewidth=1.3)
plt.xlabel('Frequency [Hz]', fontsize=25, labelpad=10)
plt.ylabel('Power', fontsize=25, labelpad=10)
plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
plt.text(0.008, 10**-0.41, r'$A$ =  {0:0.3e}'.format(m2_param[0]), fontsize=25)
plt.text(0.008, 10**-0.58, r'$n$ =  {0:0.3f}'.format(m2_param[1]), fontsize=25)
plt.text(0.008, 10**-0.75, r'$C$ =  {0:0.3e}'.format(m2_param[2]), fontsize=25)
plt.text(0.008, 10**-0.92, r'$\alpha$ =  {0:0.3f}'.format(m2_param[3]), fontsize=25)
plt.text(0.008, 10**-1.09, r'$\beta$ = {0:0.3f}'.format(m2_param[4]), fontsize=25)
plt.text(0.008, 10**-1.26, r'$\sigma$ =  {0:0.3f}'.format(m2_param[5]), fontsize=25)
plt.text(0.008, 10**-1.43, r'$\chi^2$ = {0:0.4f}'.format(redchisqrM2), fontsize=25)
plt.text(0.008, 10**-1.60, r'$r$ = {0:0.4f}'.format(r), fontsize=25)
plt.legend(loc='lower left', prop={'size':23})
#plt.show()
#plt.savefig('C:/Users/Brendan/Desktop/test_format/171_%ix_%iy_4_one.jpeg' % (m2[m],l2[m]))
#plt.savefig('C:/Users/Brendan/Desktop/171_slice2_double_optimize/171A_%ii_%ij.jpeg' % (l,m))
#plt.savefig('C:/Users/Brendan/Desktop/171_points_square/pixel_%ii_%ij_new.jpeg' % (l2[m],m2[m]))
#plt.savefig('C:/Users/Brendan/Desktop/jacks_3.0_no_averaging_index_%i_ds_std_dev.pdf' % indx, format='pdf')
#plt.close()
#"""