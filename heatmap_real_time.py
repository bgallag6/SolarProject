# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 23:35:16 2017

@author: Brendan
"""

# real-time param heatmap - 8x speed cost if show every pixel
# if show every row - 1.5x speed cost

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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sunpy
from sunpy.map import Map


#SPECTRA = np.load('C:/Users/Brendan/Desktop/SDO/spectra_20130815_193_1000_1600i_1950_2950j_rebin2.npy')
SPECTRA = np.load('C:/Users/Brendan/Desktop/SDO/spectra_20130530_1600_2300_2600i_2200_3000j_data_rebin4.npy')

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
    

## load in array of segment-averaged pixel FFTs
spectra_array = SPECTRA
#spectra_array = spectra_array.astype(np.float32)  # possibly use this?

print "The region size is %ii x %ij" % (SPECTRA.shape[0], SPECTRA.shape[1])
print "%i frequencies were evaluated in the FFT" % SPECTRA.shape[2] 

num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used
    
# determine frequency values that FFT will evaluate
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # add as argument, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]
#freqs = freqs.astype(np.float32)  # possibly use this?


# initialize arrays to hold parameter values, also each pixel's combined model fit - for tool
diffM1M2 = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1]))  # dont really use - get rid of?
#params = np.zeros((7, SPECTRA.shape[0], SPECTRA.shape[1]))

#M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], (len(freqs)+1)/2))  # would save storage / memory space
M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], SPECTRA.shape[2]))

Uncertainties = np.zeros((6, SPECTRA.shape[0], SPECTRA.shape[1]))


#20130815_193_rebin2
#x1 = 250
#x2 = 350
#y1 = 140
#y2 = 240
#bound_min = [0.1e-7, 1.5, 0.0001, 0.0005, -6.5, 0.08]
#bound_max = [2.5e-7, 2.35, 0.0015, 0.0055, -4.6, 0.8]   

#20130530_1600_rebin4
x1 = 90
x2 = 115
y1 = 45
y2 = 70
bound_min = [0.1e-7, 1.5, -0.000015, 0.008, -6.1, 0.34]
bound_max = [1.6e-7, 2.1, 0.0000025, 0.017, -5.5, 0.47]   

params = np.zeros((7, (y2-y1), (x2-x1)))  # just for real-time demo

#vis = np.load('F:/Users/Brendan/Desktop/SolarProject/visual/visual_20130815_193_1000_1600i_1950_2950j.npy')
vis = np.load('E:/Users/Brendan/Desktop/SolarProject/visual/visual_20130530_1600_2300_2600i_2200_3000j_data_rebin4.npy')
visual = vis[0]
visual = visual[y1:y2,x1:x2]
#visual = visual[x1:x2,y1:y2]

titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width', '$/chi^2$', 'Visual Image - Averaged']
   
date = '2013/05/30'
wavelength = 1600
par = 1
    
        
v_min = np.percentile(visual,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
v_max = np.percentile(visual,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)     

fig = plt.figure(figsize=(20,10))
#fig, ax = plt.subplots(1, 1, figsize=(20,15))
ax = plt.subplot2grid((1,11),(0, 0), colspan=5, rowspan=1)
ax = plt.gca()
plt.title(r'%s: %i $\AA$ -- Visual: Average' % (date, wavelength), y = 1.01, fontsize=25)
#im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
im = ax.imshow(visual, cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
#cbar.set_label('Intensity', size=20, labelpad=10)
cbar.ax.tick_params(labelsize=13, pad=3) 


start = timer()
T1 = 0



#fig, ax1 = plt.subplots(1, 1, figsize=(20,15))
ax1 = plt.subplot2grid((1,11),(0, 6), colspan=5, rowspan=1)
ax1 = plt.gca()
canvas = ax1.figure.canvas#
ax1.set_title('%s : %i $\AA$ -- %s' % (date, wavelength, titles[par]), y = 1.01, fontsize=23)
im = ax1.imshow(params[par], vmin=bound_min[par], vmax=bound_max[par])
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
cbar.ax.tick_params(labelsize=13, pad=3) 
#ax.hold(True)


plt.show(False)#
plt.draw()#
# cache the background
background = fig.canvas.copy_from_bbox(ax1.bbox)
fig.canvas.draw()


### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

#for l in range(0,SPECTRA.shape[0]):
for l in range(y1,y2):
    
    #for m in range(0,SPECTRA.shape[1]):
    for m in range(x1,x2):
        
                                        
        f = freqs  # frequencies
        s = spectra_array[l][m]  # fourier power
        
        #ds = (1./f**2.2)/1000
        ds = s*0.1  # set the error / variance estimate to a constant percentage of the spectra power-values
        
        # create points to fit model with final parameters 
        #f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],(len(freqs)+1)/2)  # would save storage / memory space?
        f_fit = freqs       
        
                                               
        ### fit data to models using SciPy's Levenberg-Marquart method
        
        try:
            # initial guesses for fitting parameters
            M1_low = [-0.002, 0.3, -0.01]
            M1_high = [0.002, 4., 0.01]
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
            M2_high = [0.002, 4., 0.01, 0.2, -4.6, 0.8]
            #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
            
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
        Uncertainties[:, l, m] = uncertainties_arr
        
        
        # create model functions from fitted parameters
        m1_fit = PowerLaw(f_fit, A, n, C)
        m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
        s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
        m2P_fit = PowerLaw(f_fit, A2, n2, C2)
        m2G_fit = Gauss(f_fit, P2, fp2, fw2)
        
        diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
        diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
                               
        
        residsM2 = (s - s_fit_gp_full)
        chisqrM2 = ((residsM2/ds)**2).sum()
        redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
        
        residsM1 = (s - m1_fit)
        chisqrM1 =  ((residsM1/ds)**2).sum()
        redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       
        
        f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))
        
        """
        # populate array with parameters
        params[0][l][m] = A2
        params[1][l][m] = n2
        params[2][l][m] = C2
        params[3][l][m] = P2
        params[4][l][m] = fp2
        params[5][l][m] = fw2
        #params[6][l][m] = redchisqrM2
        params[6][l][m] = f_test
        """
        
        params[0][l-y1][m-x1] = A2
        params[1][l-y1][m-x1] = n2
        params[2][l-y1][m-x1] = C2
        params[3][l-y1][m-x1] = P2
        params[4][l-y1][m-x1] = fp2
        params[5][l-y1][m-x1]= fw2
        #params[6][l][m] = redchisqrM2
        params[6][l-y1][m-x1] = f_test
        
        # populate array holding model fits
        M2_fit[l][m] = m2_fit
        
        ### plot only updates (pixel by pixel)
        # restore background
        canvas.restore_region(background)
        # redraw just the points
        im.set_data(params[par])
        # fill in the figure
        canvas.blit(ax1.bbox)
        plt.pause(0.0001)
        
    ### plot only updates (row by row)    
    # restore background
    #canvas.restore_region(background)
    # redraw just the points
    #im.set_data(params[par])
    # fill in the figure
    #canvas.blit(ax1.bbox)
    #plt.pause(0.0001)

    # estimate time remaining and print to screen  (looks to be much better - not sure why had above?)
    T = timer()
    T2 = T - T1
    if l == 0:
        T_init = T - start
        T_est = T_init*(SPECTRA.shape[0])  
    else:
        T_est2 = T2*((SPECTRA.shape[0])-l)
        print "Currently on row %i of %i, estimated time remaining: %i seconds" % (l, SPECTRA.shape[0], T_est2)
    T1 = T

# print estimated and total program time to screen        
T_act = timer() - start
print "Actual total time = %i sec" % T_act           

