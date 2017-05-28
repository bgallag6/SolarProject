# -*- coding: utf-8 -*-
"""
Created on Thu May 18 12:38:39 2017

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
from scipy.stats import f as ff
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
    
directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20130815'
wavelength = 171

spec_std = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_std.npy')
spec_segs = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_9_test.npy')

#spec1 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/spec_array1_jack.npy')
#spec3 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/spec_array3_jack.npy')
#spec6 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/spec_array6_jack.npy')
#spec12 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/spec_array12_jack.npy')

"""
spec1 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/spec_array_1.npy')
spec3 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/spec_array_3.npy')
spec6 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/spec_array_6.npy')
spec12 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/spec_array_12.npy')

stddev1 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/std_dev_1.npy')
stddev3 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/std_dev_3.npy')
stddev6 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/std_dev_6.npy')
stddev12 = np.load('C:/Users/Brendan/Desktop/solar/spectra_fitting/spec_normalized/std_dev_12.npy')
"""

spec1 = np.load('C:/Users/Brendan/Desktop/spec_array_1.npy')
spec3 = np.load('C:/Users/Brendan/Desktop/spec_array_3.npy')
spec6 = np.load('C:/Users/Brendan/Desktop/spec_array_6.npy')
spec12 = np.load('C:/Users/Brendan/Desktop/spec_array_12.npy')

stddev1 = np.load('C:/Users/Brendan/Desktop/std_dev_1.npy')
stddev3 = np.load('C:/Users/Brendan/Desktop/std_dev_3.npy')
stddev6 = np.load('C:/Users/Brendan/Desktop/std_dev_6.npy')
stddev12 = np.load('C:/Users/Brendan/Desktop/std_dev_12.npy')
    
#spectra_array = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/193/20130626_193_-450_-200i_-200_200j_spectra.npy')
#spectra_array = np.load('C:/Users/Brendan/Desktop/project_files/20130626_171_-500_500i_-500_600j_spectra_arth.npy')
cube_shape = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength))
spectra_array = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))
#spectra_array = np.load('C:/Users/Brendan/Desktop/20160905/193_1000.npy')
#spectra_array = np.load('C:/Users/Brendan/Desktop/1600/spectra.npy')
## load in array of segment-averaged pixel FFTs
SPECTRA = spectra_array

## determine frequency values that FFT will evaluate
if wavelength == 1600 or wavelength == 1700:
  #time_step = 24  # 24-second cadence for these wavelengths
  time_step = 12  # Jacks
else:
  time_step = 12  # 12-second cadence for the others
  
num_seg = np.array([1,3,6,12])   
num_freq = np.array([1799,599,299,149])
specs = np.array([spec1,spec3,spec6,spec12]) 
stddev = np.array([stddev1,stddev3,stddev6,stddev12])   

A0 = 6e-9
n0 = 1.6
C0 = 2e-04

for n in range(len(num_seg)):
#for n in range(1):
 
    n_segments = num_seg[n]  # break data into 12 segments of equal length
    n_freqs = num_freq[n]
    
    freq_size = ((n_freqs)*2) + 1  # determined from FFT-averaging script
    time_step = 12  # add as argument, or leave in as constant?
    #time_step = 24  # add as argument, or leave in as constant?
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    if n_segments == 1:
        freqs1 = freqs
    elif n_segments == 3:
        freqs3 = freqs  
    elif n_segments == 6:
        freqs6 = freqs
    elif n_segments == 12:
        freqs12 = freqs

    f = freqs
    s = specs[n]
    
    # assign equal weights to all parts of the curve
    #df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
    #df2 = np.zeros_like(f)
    #df2[0:len(df)] = df
    #df2[len(df2)-1] = df2[len(df2)-2]
    #ds = df2
    #ds = 0.1*s
    #ds = spec_std
    ds = stddev[n]    
    #plt.loglog(f,ds)
    #plt.loglog(f,s)
    
    M1_low = [-0.002, 0.3, -0.01]
    M1_high = [0.002, 6., 0.01]
    nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
       
    
     
    A, n, C = nlfit_l  # unpack fitting parameters
    
    # unpack uncertainties in fitting parameters from diagonal of covariance matrix
    dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
    
    #gauss_wid = 8                            
    M2_low = [-0.002, 0.3, -0.01, 0.000001, -6.5, 0.05]
    #M2_low = [-0.002, 0.3, -0.01, -0.2, -6.5, 0.05]
    M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
    #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
    
    # change method to 'dogbox' and increase max number of function evaluations to 3000
    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.43], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
        
    
    
    A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
    
    # unpack uncertainties in fitting parameters from diagonal of covariance matrix
    dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
    
        
    
    # create model functions from fitted parameters
    #m1_fit = PowerLaw(f_fit, A, n, C)
    m1_fit = PowerLaw(f, A, n, C)
    amp_scale = PowerLaw(np.exp(fp2), A, n, C)  # to extract the gaussian-amplitude scaling factor
    #m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
    m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
    #s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
    
    """
    m2P_fit = PowerLaw(f, A2, n2, C2)  # only need if plotting
    m2G_fit = Gauss(f, P2, fp2, fw2)  # only need if plotting
    m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]      
    uncertainties = dA2, dn2, dC2, dP2, dfp2, dfw2  # do we want to keep a global array of uncertainties?
    """
    
    #diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
    #diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
    #M2_low = [0., 0.3, 0., 0., -6.5, 0.05]
    #M2_high = [0.002, 6., 0.01, 0.1, -4.6, gauss_wid/10.]
    #"""
    nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
    #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
    A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
    dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
    m2_param = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
    m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22) 
    uncertainties = dA22, dn22, dC22, dP22, dfp22, dfw22  # do we want to keep a global array of uncertainties?
    
    m2P_fit = PowerLaw(f, A22, n22, C22)  # only need if plotting
    m2G_fit = Gauss(f, P22, fp22, fw22)  # only need if plotting    
    
    r_val = pearsonr(m2_fit2, s)
    
    if n_segments == 1:
        m2_fit1 = m2_fit2
        m2_param1 = A22, n22, C22, P22, fp22, fw22
        m2G_fit1 = m2G_fit
        r_val1 = r_val[0]
    elif n_segments == 3:
        m2_fit3 = m2_fit2  
        m2_param3 = A22, n22, C22, P22, fp22, fw22
        m2G_fit3 = m2G_fit
        r_val3 = r_val[0]
    elif n_segments == 6:
        m2_fit6 = m2_fit2
        m2_param6 = A22, n22, C22, P22, fp22, fw22
        m2G_fit6 = m2G_fit 
        r_val6 = r_val[0]
    elif n_segments == 12:
        m2_fit12 = m2_fit2
        m2_param12 = A22, n22, C22, P22, fp22, fw22
        m2G_fit12 = m2G_fit 
        r_val12 = r_val[0]
    
    residsM22 = (s - m2_fit2)
    chisqrM22 = ((residsM22/ds)**2).sum()
    redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6)  
    
           
    #"""        
    
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
        
    #f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))
    f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))
    
    #amp_scale = PowerLaw(np.exp(fp2), A2, n2, C2)  # to extract the gaussian-amplitude scaling factor
    amp_scale2 = PowerLaw(np.exp(fp22), A22, n22, C22)  # to extract the gaussian-amplitude scaling factor
    
    # generate p-value heatmap
    df1, df2 = 3, 6
    p_val = ff.sf(f_test2, df1, df2)
    
    
    
    """
    fig = plt.figure(figsize=(15,15))
    ax = plt.gca()  # get current axis -- to set colorbar 
    #plt.title('Power-Law Dominated : Pixel %ii, %ij' % (l2[m],m2[m]), y = 1.01, fontsize=25)
    plt.title('%s: Pixel %ix, %iy' % (m2_title[m], m2[m],l2[m]), y = 1.01, fontsize=30)
    plt.ylim((10**-4.7,10**0))
    plt.xlim((10**-4.,10**-1.3))
    plt.xticks(fontsize=19)
    plt.yticks(fontsize=19)
    plt.loglog(f,s,'k')
    plt.loglog(f, m1_fit, label='M1 - Power Law', linewidth=1.3)
    plt.loglog(f, m2P_fit, 'g', label='M2 - Power Law', linewidth=1.3)
    plt.loglog(f, m2G_fit, 'g--', label='M2 - Gaussian', linewidth=1.3)
    plt.loglog(f, m2_fit2, 'purple', label='M2 - Combined', linewidth=1.3)
    plt.xlabel('Frequency [Hz]', fontsize=25, labelpad=10)
    plt.ylabel('Power', fontsize=25, labelpad=10)
    plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
    plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
    
    #rect = patches.Rectangle((0.005,0.05), 0.03, 0.6, color='white', fill=True)
    #ax.add_patch(rect)
    plt.text(0.008, 10**-0.38, r'$A$ =  {0:0.3e}'.format(m2_param[0]), fontsize=25)
    plt.text(0.008, 10**-0.58, r'$n$ =  {0:0.3f}'.format(m2_param[1]), fontsize=25)
    #plt.text(0.008, 10**-0.75, r'$C$ =  {0:0.3e}'.format(m2_param[2]), fontsize=25)
    plt.text(0.008, 10**-0.75, r'$(C/A)^{-\frac{1}{n}}$ = %.3e Seconds' % m2_param[2], fontsize=25)
    plt.text(0.008, 10**-0.92, r'$\alpha$ =  {0:0.3f}'.format(m2_param[3]), fontsize=25)
    plt.text(0.008, 10**-1.09, r'$\beta$ = {0:0.3f}'.format(m2_param[4]), fontsize=25)
    plt.text(0.008, 10**-1.26, r'$\sigma$ =  {0:0.3f}'.format(m2_param[5]), fontsize=25)
    plt.text(0.008, 10**-1.43, r'$\chi^2$ = {0:0.4f}'.format(redchisqrM22), fontsize=25)
    plt.text(0.008, 10**-1.60, r'$P-Value$ = {0:0.4f}'.format(p_val), fontsize=25)
    plt.legend(loc='lower left', prop={'size':25})
    #plt.show()
    plt.savefig('C:/Users/Brendan/Desktop/test_format2/171_%ix_%iy_font_25.pdf' % (m2[m],l2[m]), format='pdf')
    #plt.savefig('C:/Users/Brendan/Desktop/171_slice2_double_optimize/171A_%ii_%ij.jpeg' % (l,m))
    #plt.savefig('C:/Users/Brendan/Desktop/171_points_square/pixel_%ii_%ij_new.jpeg' % (l2[m],m2[m]))
    #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
    #plt.close()
    """
plt.rcParams["font.family"] = "Times New Roman"
font_size = 27

fig = plt.figure(figsize=(15,15))
ax = plt.gca()  # get current axis -- to set colorbar 
#plt.title('Comparison of Temporal Averaging Methods - Normalized by # of Segments \n Synthetic Parameters : A = %0.2g | n = %0.2f | R = %0.1f | ds = log10(f+1) - log10(f)' % (A0,n0,((C0/A0)**(1./n0))/60), y = 1.01, fontsize=25)
plt.title('Comparison of Temporal Averaging Methods - Normalized by # of Segments', y = 1.01, fontsize=25)
#plt.title('Comparison of Temporal Averaging Methods \n Unnormalized by # of Segments', y=1.01, fontsize=23)
#plt.title('%s: Pixel %ix, %iy' % (m2_title[m], m2[m],l2[m]), y = 1.01, fontsize=30)
#plt.title('171: %ix, %iy' % (m2[m],l2[m]), y = 1.01, fontsize=30)
#plt.title('%s: Point %s' % (m2_title[m], point_label[m]), y = 1.01, fontsize=font_size, fontname="Times New Roman")
#plt.title('%s: Point %s' % (m2_title[m-17], point_label[m-17]), y = 1.01, fontsize=30)
#plt.ylim((10**-4.7,10**0))
#plt.xlim((10**-4.,10**-1.3))
plt.ylim(10**-6.5,10**0)
plt.xlim(10**-5.,10**-1.3)
plt.xticks(fontsize=font_size, fontname="Times New Roman")
plt.yticks(fontsize=font_size, fontname="Times New Roman")
ax.tick_params(axis='both', which='major', pad=10)

#plt.loglog(freqs,spec_segs[4])
#plt.loglog(f,s,'k',linewidth=2, label='Synthetic + Noise')
#plt.loglog(f, temp_avg)
#plt.loglog(f, temp_avg2)
#plt.loglog(f,ds,'r',linewidth=2)
#plt.loglog(f,spec_std,'purple',linewidth=2)
#plt.loglog(f, m1_fit, label='M1 - Power Law', linewidth=1.3)
#plt.loglog(f, m2P_fit, 'g', label='M2 - Power Law', linewidth=1.3)
#plt.loglog(f, m2G_fit, 'g--', label='M2 - Gaussian', linewidth=1.3)
plt.loglog(freqs1, m2_fit1, 'k', linewidth=1.3, label='  (1) 12-Hour Segment: A = %0.3g | n = %0.2f | R = %0.1f   | r = %0.3f' % (m2_param1[0],m2_param1[1],((m2_param1[2]/m2_param1[0])**(1./m2_param1[1]))/60, r_val1))
#plt.loglog(freqs3, m2_fit3, 'b', linewidth=1.3, label='  (3) 4-Hour Segments: A = %0.3g  | n = %0.2f | R = %0.1f   | r = %0.3f' % (m2_param3[0],m2_param3[1],((m2_param3[2]/m2_param3[0])**(1./m2_param3[1]))/60, r_val3))
plt.loglog(freqs6, m2_fit6, 'g', linewidth=1.3, label='  (6) 2-Hour Segments: A = %0.3g | n = %0.2f | R = %0.1f | r = %0.3f' % (m2_param6[0],m2_param6[1],((m2_param6[2]/m2_param6[0])**(1./m2_param6[1]))/60, r_val6))
#plt.loglog(freqs12, m2_fit12, 'r', linewidth=1.3, label='(12) 1-Hour Segments: A = %0.3g | n = %0.2f | R = %0.1f | r = %0.3f' % (m2_param12[0],m2_param12[1],((m2_param12[2]/m2_param12[0])**(1./m2_param12[1]))/60, r_val12))
plt.loglog(freqs1, spec1, 'k', linewidth=1.3)
#plt.loglog(freqs3, spec3, 'b', linewidth=1.3)
plt.loglog(freqs6, spec6, 'g', linewidth=1.3)
#plt.loglog(freqs12, spec12, 'r', linewidth=1.3)
plt.loglog(freqs1, m2G_fit1, 'k', linewidth=1.3)
#plt.loglog(freqs3, m2G_fit3, 'b', linewidth=1.3)
plt.loglog(freqs6, m2G_fit6, 'g', linewidth=1.3)
#plt.loglog(freqs12, m2G_fit12, 'r', linewidth=1.3)
#plt.loglog(f, m1_fit, linewidth=1.3)
#plt.loglog(f, m2P_fit, 'g', label=r'$A\nu^{-n} + C$', linewidth=1.3)
#plt.loglog(f, m2G_fit, 'g--', label=r'$\alpha\ e^{{-\frac{(\ln\nu-\beta)^{2}}{\sigma^{2}}}}$', linewidth=1.3)
#plt.loglog(f, m2_fit2, 'purple', label='Combined Model', linewidth=1.3)
plt.xlabel('Frequency [Hz]', fontsize=font_size, labelpad=10, fontname="Times New Roman")
plt.ylabel('Power', fontsize=font_size, labelpad=10, fontname="Times New Roman")
plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed')
plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted')
#plt.vlines((0.0093),10**-8,10**1, linestyles='dotted', label='3 minutes')

#rect = patches.Rectangle((0.004,0.015), 0.012, 0.8, color='white', fill=True)
#ax.add_patch(rect)
#"""
#plt.text(0.0036, 10**-0.91, r'Gaussian Width = %0.2f' % (gauss_wid/10.), fontsize=font_size, fontname="Times New Roman")
#plt.text(0.0076, 10**-0.61, r'$n$ = {0:0.2f}'.format(m2_param[1]), fontsize=font_size, fontname="Times New Roman")
#plt.text(0.008, 10**-0.75, r'$C$ =  {0:0.3e}'.format(m2_param[2]), fontsize=25)
#plt.text(0.007, 10**-0.73, r'$(C/A)^{-\frac{1}{n}}$ = %f [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
#plt.text(0.00725, 10**-0.83, r'$R$ = %0.1f [min]' % ((1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1]))/60.), fontsize=font_size, fontname="Times New Roman")
#plt.text(0.007, 10**-0.73, r'$r$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
#plt.text(0.0075, 10**-1.05, r'$\alpha$ = {0:0.2e}'.format(m2_param[3]), fontsize=font_size, fontname="Times New Roman")
#plt.text(0.007, 10**-1.09, r'$\beta$ = {0:0.3f}'.format(m2_param[4]), fontsize=25)
#plt.text(0.00763, 10**-1.26, r'$\beta$ = {0:0.1f} [min]'.format((1./np.exp(m2_param[4]))/60.), fontsize=font_size, fontname="Times New Roman")
#plt.text(0.0079, 10**-1.47, r'$\delta$ = {0:0.3f}'.format(m2_param[5]), fontsize=font_size, fontname="Times New Roman")
#plt.text(0.007, 10**-1.55, r'$\chi^2$ = {0:0.3f}'.format(chisqrM22), fontsize=30)
#plt.text(0.007, 10**-1.75, r'$p$ = {0:0.2e}'.format(p_val), fontsize=30)
#plt.text(0.0077, 10**-1.67, r'$p$ = {0:0.3g}'.format(p_val), fontsize=font_size, fontname="Times New Roman")
#plt.text(0.00786, 10**-1.87, r'$r$ = {0:0.3g}'.format(r_val[0]), fontsize=font_size, fontname="Times New Roman")
#legend = ax.legend(loc='upper right', prop={'size':30}, labelspacing=0.35)
legend = ax.legend(loc='lower left', prop={'size':25}, labelspacing=0.35)
for label in legend.get_lines():
    label.set_linewidth(3.0)  # the legend line width
#"""
"""
plt.text(0.00015, 10**-3.15, r'$A$ =  {0:0.3e}'.format(m2_param[0]), fontsize=25)
plt.text(0.00015, 10**-3.32, r'$n$ =  {0:0.3f}'.format(m2_param[1]), fontsize=25)
#plt.text(0.008, 10**-0.75, r'$C$ =  {0:0.3e}'.format(m2_param[2]), fontsize=25)
plt.text(0.00015, 10**-3.52, r'$(C/A)^{-\frac{1}{n}}$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=25)
plt.text(0.00015, 10**-3.72, r'$\alpha$ =  {0:0.3f}'.format(m2_param[3]), fontsize=25)
plt.text(0.00015, 10**-3.89, r'$\beta$ = {0:0.3f}'.format(m2_param[4]), fontsize=25)
plt.text(0.00015, 10**-4.06, r'$\sigma$ =  {0:0.3f}'.format(m2_param[5]), fontsize=25)
plt.text(0.00015, 10**-4.23, r'$\chi^2$ = {0:0.4f}'.format(redchisqrM22), fontsize=25)
plt.text(0.00015, 10**-4.42, r'$P-Value$ = {0:0.3f}'.format(p_val), fontsize=25)
plt.legend(loc='upper right', prop={'size':23})
"""
letter = 'C'
#plt.savefig('C:/Users/Brendan/Desktop/spec_avg_fits_arithmetic%s.pdf' % letter, format='pdf')
#plt.savefig('C:/Users/Brendan/Desktop/spec_avg_fits_geometric%s.pdf' % letter, format='pdf')
