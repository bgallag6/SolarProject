# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 14:35:48 2017

@author: Brendan
"""

import numpy as np
import scipy.signal
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
from scipy import fftpack  # doesnt work in module when called here???
from astropy.convolution import convolve, Box1DKernel
from numpy.random import randn
from mpi4py import MPI
import matplotlib.pyplot as plt

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



spectra = np.load('C:/Users/Brendan/Desktop/validation/1seg_1x1_spectra.npy')   


num_freq = spectra.shape[1]  # determine nubmer of frequencies that are used
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # 12 second cadence for the others
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

a0 = 10**-9.
n0 = 2.
c0 = 10**-4.5
p0 = 10**-3.5
l0 = -5.5
w0 = 0.3

ex = 'L'

"""
# maybe put how the fit changes depending on the parameters in the slider tool?
"""

spectra = a0*freqs**(-n0)+c0+p0*np.exp(-0.5*(((np.log(freqs))-(l0))/w0)**2)

#noise1 = np.random.ranf(len(spectra))
noise1 = np.random.uniform(-1.0,1.0,len(spectra))
#plt.figure()
#plt.hist(noise1)

sig = spectra

#spectra = sig+(sig*0.5*noise1)
spectra = sig+(sig*0.5*noise1)

orig_powerspec = spectra

# Now Smooth the Power Spectra
orig = orig_powerspec.copy()		# make a copy

window_len=11   # Keep this an odd number

window='hanning'  # OTHER CHOICES =  'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
s=np.r_[orig[window_len-1:0:-1],orig,orig[-1:-window_len:-1]]

if window == 'flat': #moving average
	w=np.ones(window_len,'d')
else:
	w=eval('np.'+window+'(window_len)')

	
y=np.convolve(w/w.sum(),s,mode='valid')
y=y[( window_len/2+1): len(orig)+( window_len/2 +1)]	# Crop down the new TS

# Now smooth the smooth...
y2c = y.copy()
s2=np.r_[y2c[window_len-1:0:-1],y2c,y2c[-1:-window_len:-1]]
y2=np.convolve(w/w.sum(),s2,mode='valid')
y2=y2[( window_len/2+1): len(orig)+( window_len/2 +1)]   # Crop down the new TS


# PLOTS TO SEE HOW IT WORKED OUT
plt.figure(figsize=(15,15))
plt.title('Example %s \n Hanning-Window Smoothing' % (ex), y = 1.01, fontsize=30)
plt.loglog(freqs, sig,color='blue', label='Pure Signal')
plt.loglog(freqs, orig,color='black', label='Spectra w/ Noise')
plt.loglog(freqs, y,color='green', label='1st Hanning Smooth')
plt.loglog(freqs, y2,color='red', label='2nd Hanning Smooth')
ax = plt.gca()  # get current axis -- to set colorbar 
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlim(10**-4.5,10**-1)
plt.ylim(10**-5.5,10**-.01)
ax.tick_params(axis='both', which='major', pad=15)
plt.xlabel('Frequency [Hz]', fontsize=30, labelpad=10)
plt.ylabel('Power', fontsize=30, labelpad=10)
plt.legend(loc='upper right')
#plt.savefig('C:/Users/Brendan/Desktop/param_fits/example_hanning_%s.pdf' % ex, format='pdf')

"""
# determine frequency values that FFT will evaluate
num_freq = spectra.shape[1]  # determine nubmer of frequencies that are used
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # 12 second cadence for the others
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

# initialize arrays to hold parameter values, also each pixel's combined model fit - for tool
# M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], (len(freqs)+1)/2))  # would save storage / memory space
#M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], SPECTRA.shape[2]))     
params = np.zeros((8, spectra.shape[0], spectra.shape[1]))
"""
#for ii in range(spectra.shape[0]):  
for ii in range(1):  
    print ii     
    #for jj in range(spectra.shape[0]): 
    for jj in range(1):                             
        f = freqs  # frequencies
        #s = spectra[ii][jj]
        #s = spectra[ii]
        s = y2
        
        #plt.figure()
        #plt.loglog(f,s)
        
        
        # assign equal weights to all parts of the curve
        df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
        df2 = np.zeros_like(f)
        df2[0:len(df)] = df
        df2[len(df2)-1] = df2[len(df2)-2]
        ds = df2
        
        # create points to fit model with final parameters 
        #f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],(len(freqs)+1)/2)  # would save storage / memory space?      
        
                                               
        ### fit data to models using SciPy's Levenberg-Marquart method
        
        try:
            # initial guesses for fitting parameters
            M1_low = [-0.002, 0.3, -0.01]
            M1_high = [0.002, 6., 0.01]
            nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
           
        
        except RuntimeError:
            #print("Error M1 - curve_fit failed - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        except ValueError:
            #print("Error M1 - inf/NaN - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
          
        A, n, C = nlfit_l  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
        
        ## fit data to combined power law plus gaussian component model
               
        try:                                 
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
            #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
            
            # change method to 'dogbox' and increase max number of function evaluations to 3000
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='trf', max_nfev=3000) # replaced #'s with arrays
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.425], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
        
        except RuntimeError:
            #print("Error M2 - curve_fit failed - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        except ValueError:
            #print("Error M2 - inf/NaN - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
        
        try:
            nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
            #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
           
        except RuntimeError:
            #print("Error M2 - curve_fit failed - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        except ValueError:
            #print("Error M2 - inf/NaN - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
        dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
                       
        # create model functions from fitted parameters
        m1_fit = PowerLaw(f, A, n, C)        
        #m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
        m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22)   
        m2P_fit = PowerLaw(f, A22, n22, C22)  # only need if plotting
        m2G_fit = Gauss(f, P22, fp22, fw22)  # only need if plotting   
        m2_param = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
        
        residsM1 = (s - m1_fit)
        chisqrM1 =  ((residsM1/ds)**2).sum()
        redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)  
                       
        #residsM2 = (s - m2_fit)
        #chisqrM2 = ((residsM2/ds)**2).sum()
        #redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
        
        residsM22 = (s - m2_fit2)
        chisqrM22 = ((residsM22/ds)**2).sum()
        redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6) 
        
        plt.figure(figsize=(15,15))
        plt.title(r"Example %s" % (ex) + "\n" + r"A: %0.2e | n: %0.2f | $R$: %1.0f | $\alpha$: %0.2e | $\beta$: %1.0f | $\sigma$: %0.3f" % (a0,n0,((c0/a0)**(1./n0)),p0,(1./np.exp(l0)),w0), y = 1.01, fontsize=30)
        plt.loglog(f,s)
        plt.loglog(f,m2_fit2)
        ax = plt.gca()  # get current axis -- to set colorbar 
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlim(10**-4.5,10**-1)
        plt.ylim(10**-5.5,10**-0.01)
        ax.tick_params(axis='both', which='major', pad=15)
        plt.xlabel('Frequency [Hz]', fontsize=30, labelpad=10)
        plt.ylabel('Power', fontsize=30, labelpad=10)
        
        #rect = patches.Rectangle((0.005,0.05), 0.03, 0.6, color='white', fill=True)
        #ax.add_patch(rect)
        
        plt.text(0.0045, 10**-0.51, r'$A$ = {0:0.2e}'.format(m2_param[0]), fontsize=30)
        plt.text(0.0046, 10**-0.81, r'$n$ = {0:0.2f}'.format(m2_param[1]), fontsize=30)
        #plt.text(0.008, 10**-0.75, r'$C$ =  {0:0.3e}'.format(m2_param[2]), fontsize=25)
        #plt.text(0.007, 10**-0.73, r'$(C/A)^{-\frac{1}{n}}$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
        plt.text(0.0045, 10**-1.13, r'$R$ = %1.0f [$s$]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
        #plt.text(0.007, 10**-0.73, r'$r$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
        plt.text(0.0046, 10**-1.45, r'$\alpha$ = {0:0.2e}'.format(m2_param[3]), fontsize=30)
        #plt.text(0.007, 10**-1.09, r'$\beta$ = {0:0.3f}'.format(m2_param[4]), fontsize=25)
        plt.text(0.0046, 10**-1.75, r'$\beta$ = {0:1.0f} [$s$]'.format(1./np.exp(m2_param[4])), fontsize=30)
        plt.text(0.0046, 10**-2.05, r'$\sigma$ = {0:0.3f}'.format(m2_param[5]), fontsize=30)
        #plt.text(0.007, 10**-1.55, r'$\chi^2$ = {0:0.3f}'.format(chisqrM22), fontsize=30)
        #plt.text(0.007, 10**-1.75, r'$p$ = {0:0.2e}'.format(p_val), fontsize=30)
        #plt.text(0.0047, 10**-1.55, r'$p$ = {0:0.3g}'.format(p_val), fontsize=30)
        #plt.text(0.0047, 10**-1.75, r'$r$ = {0:0.3g}'.format(r_val[0]), fontsize=30)
        #plt.savefig('C:/Users/Brendan/Desktop/param_fits/example_hanning_fit_%s.pdf' % ex, format='pdf')
        
        """
        plt.figure(figsize=(15,15))
        plt.hist(residsM22,bins=100,range=(-0.002,0.002))
        plt.title('Residuals from 1-Segment No Pixel-Box Averaging', y = 1.01, fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel('Residual Amount', fontsize=25, labelpad=10)
        plt.ylabel('Bin Count', fontsize=25, labelpad=10)
        #plt.savefig('C:/Users/Brendan/Desktop/validation/histogram_1_1x1.pdf', format='pdf')
        """        
        """
        ds2 = np.abs(residsM22)
    
        try:
            #M2_low = [-1., 0.3, -1., 0., -6.5, 0.05]
            #M2_high = [1., 6., 1., 0.5, -4.6, 0.8]
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
            nlfit_gp3, nlpcov_gp3 = scipy.optimize.curve_fit(GaussPowerBase, f, s,  bounds=(M2_low, M2_high), sigma=ds2, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
           
        except RuntimeError:
            #print("Error M2 - curve_fit failed - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        except ValueError:
            #print("Error M2 - inf/NaN - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        
        A23, n23, C23, P23, fp23, fw23 = nlfit_gp3  # unpack fitting parameters     
        m3_param = A23, n23, C23, P23, fp23, fw23  # could have used this for params array : = params[0:6,l-1,m-1]

        m2_fit3 = GaussPowerBase(f, A23,n23,C23,P23,fp23,fw23)           
        
        
        plt.figure(figsize=(15,15))
        ax = plt.gca()  # get current axis -- to set colorbar 
        plt.title(r'1-Segment No Pixel-Box Averaging: ds = |resids|', y = 1.01, fontsize=30)
        plt.loglog(f,s)
        plt.loglog(f,m2_fit3)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.xlim(10**-4.5,10**-1)
        plt.ylim(10**-7,10**0)
        ax.tick_params(axis='both', which='major', pad=15)
        plt.xlabel('Frequency [Hz]', fontsize=30, labelpad=10)
        plt.ylabel('Power', fontsize=30, labelpad=10)
        
        plt.text(0.0045, 10**-0.51, r'$A$ = {0:0.2e}'.format(m3_param[0]), fontsize=30)
        plt.text(0.0046, 10**-0.81, r'$n$ = {0:0.2f}'.format(m3_param[1]), fontsize=30)
        #plt.text(0.008, 10**-0.75, r'$C$ =  {0:0.3e}'.format(m2_param[2]), fontsize=25)
        #plt.text(0.007, 10**-0.73, r'$(C/A)^{-\frac{1}{n}}$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
        plt.text(0.0045, 10**-1.13, r'$R$ = %1.0f [$s$]' % (1./(m3_param[2] / m3_param[0])**(-1./ m3_param[1])), fontsize=30)
        #plt.text(0.007, 10**-0.73, r'$r$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
        plt.text(0.0046, 10**-1.45, r'$\alpha$ = {0:0.2e}'.format(m3_param[3]), fontsize=30)
        #plt.text(0.007, 10**-1.09, r'$\beta$ = {0:0.3f}'.format(m2_param[4]), fontsize=25)
        plt.text(0.0046, 10**-1.75, r'$\beta$ = {0:1.0f} [$s$]'.format(1./np.exp(m3_param[4])), fontsize=30)
        plt.text(0.0046, 10**-2.05, r'$\sigma$ = {0:0.3f}'.format(m3_param[5]), fontsize=30)
        #plt.text(0.007, 10**-1.55, r'$\chi^2$ = {0:0.3f}'.format(chisqrM22), fontsize=30)
        #plt.text(0.007, 10**-1.75, r'$p$ = {0:0.2e}'.format(p_val), fontsize=30)
        #plt.text(0.0047, 10**-1.55, r'$p$ = {0:0.3g}'.format(p_val), fontsize=30)
        #plt.text(0.0047, 10**-1.75, r'$r$ = {0:0.3g}'.format(r_val[0]), fontsize=30)
        #plt.title('1-Segment + 1x1 Pixel Box Averaging', y=1.01, fontsize=20)
        plt.savefig('C:/Users/Brendan/Desktop/validation2/fit_1_1x1_resids.pdf', format='pdf')
        
        plt.figure(figsize=(15,15))
        ax = plt.gca()  # get current axis -- to set colorbar 
        plt.title('|Residuals|', y = 1.01, fontsize=30)
        plt.loglog(f,ds2)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.xlim(10**-4.5,10**-1)
        plt.ylim(10**-7,10**0)
        ax.tick_params(axis='both', which='major', pad=15)
        plt.xlabel('Frequency [Hz]', fontsize=30, labelpad=10)
        plt.ylabel('Power', fontsize=30, labelpad=10)
        plt.savefig('C:/Users/Brendan/Desktop/validation2/resids.pdf', format='pdf')
        """
