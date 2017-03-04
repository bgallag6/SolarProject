# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 14:14:33 2017

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
spectra_array = np.load('C:/Users/Brendan/Desktop/project_files/20130626_171_-500_500i_-500_600j_spectra_arth.npy')
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

#m2 = [742, 790, 586, 572, 475, 482, 839, 842, 767, 797, 295, 134, 43, 116, 746, 746, 757, 765, 867, 867, 606, 453, 1463, 804, 828, 55, 742, 573, 525, 664, 762, 649, 711, 722, 1473, 821, 826, 1325, 408, 485, 483, 647, 743, 708, 744]
#l2 = [322, 235, 311, 313, 297, 223, 489, 590, 547, 626, 758, 371, 376, 196, 331, 325, 319, 325, 864, 807, 1124, 51, 1403, 653, 1163, 321, 323, 157, 551, 330, 1529, 1548, 1440, 1441, 1401, 1229, 1166, 1238, 427, 535, 225, 212, 322, 352, 272]

#m2 = [790, 767, 757, 765, 867, 525, 762, 649, 722, 485, 743, 708, 744, 14] #use
#l2 = [235, 547, 319, 325, 864, 551, 1529, 1548, 1441, 535, 322, 352, 272, 330] #use

#m2 = [722]
#l2 = [1441]

#m2 = [14, 58, 283, 512, 629, 810, 909, 901, 1342, 767, 873]
#l2 = [330, 394, 482, 538, 379, 293, 400, 409, 1239, 1160, 1261]

#m2 = [743, 708, 525, 757, 765, 722, 867]
#l2 = [322, 352, 551, 319, 325, 1441, 864]

m2 = [722, 525, 757, 743]
l2 = [1441, 551, 319, 322]

point_label = ['A', 'B', 'C', 'D']

m2_title = ['Tail Dominated w/o Gaussian', 'Power-Law Dominated w/o Gaussian', 'Power-Law Dominated w/ Gaussian', 'Tail Dominated w/ Gaussian']



#for l in range(321,322):
for l in range(1):
    
    #for m in range(500,900):
    #for m in range(900,901):
    for m in range(len(m2)):    
        
                                        
        f = freqs  # frequencies
        #s = spectra_array[l][m]  # fourier power
        s = spectra_array[l2[m]][m2[m]]  # fourier power
        
       
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
        
        #"""
        nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
        #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
        dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        m2_param = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
        m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22) 
        uncertainties = dA22, dn22, dC22, dP22, dfp22, dfw22  # do we want to keep a global array of uncertainties?
        
        residsM22 = (s - m2_fit2)
        chisqrM22 = ((residsM22/ds)**2).sum()
        redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6)  
        
        m2P_fit = PowerLaw(f, A22, n22, C22)  # only need if plotting
        m2G_fit = Gauss(f, P22, fp22, fw22)  # only need if plotting           
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
        
        fig = plt.figure(figsize=(15,15))
        ax = plt.gca()  # get current axis -- to set colorbar 
        #plt.title('Power-Law Dominated : Pixel %ii, %ij' % (l2[m],m2[m]), y = 1.01, fontsize=25)
        #plt.title('%s: Pixel %ix, %iy' % (m2_title[m], m2[m],l2[m]), y = 1.01, fontsize=30)
        plt.title('%s: Point %s' % (m2_title[m], point_label[m]), y = 1.01, fontsize=30)
        plt.ylim((10**-4.7,10**0))
        plt.xlim((10**-4.,10**-1.3))
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        ax.tick_params(axis='both', which='major', pad=15)
        plt.loglog(f,s,'k')
        plt.loglog(f, m1_fit, label='M1 - Power Law', linewidth=1.3)
        plt.loglog(f, m2P_fit, 'g', label='M2 - Power Law', linewidth=1.3)
        plt.loglog(f, m2G_fit, 'g--', label='M2 - Gaussian', linewidth=1.3)
        plt.loglog(f, m2_fit2, 'purple', label='M2 - Combined', linewidth=1.3)
        plt.xlabel('Frequency [Hz]', fontsize=30, labelpad=10)
        plt.ylabel('Power', fontsize=30, labelpad=10)
        plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
        plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
        
        #rect = patches.Rectangle((0.005,0.05), 0.03, 0.6, color='white', fill=True)
        #ax.add_patch(rect)
        #"""
        plt.text(0.007, 10**-0.31, r'$A$ =  {0:0.2e}'.format(m2_param[0]), fontsize=30)
        plt.text(0.007, 10**-0.51, r'$n$ =  {0:0.2f}'.format(m2_param[1]), fontsize=30)
        #plt.text(0.008, 10**-0.75, r'$C$ =  {0:0.3e}'.format(m2_param[2]), fontsize=25)
        plt.text(0.007, 10**-0.73, r'$(C/A)^{-\frac{1}{n}}$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
        #plt.text(0.007, 10**-0.73, r'$r$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
        plt.text(0.007, 10**-0.95, r'$\alpha$ =  {0:0.2e}'.format(m2_param[3]), fontsize=30)
        #plt.text(0.007, 10**-1.09, r'$\beta$ = {0:0.3f}'.format(m2_param[4]), fontsize=25)
        plt.text(0.007, 10**-1.15, r'$\beta$ = {0:1.0f} [s]'.format(1./np.exp(m2_param[4])), fontsize=30)
        plt.text(0.007, 10**-1.35, r'$\sigma$ =  {0:0.3f}'.format(m2_param[5]), fontsize=30)
        plt.text(0.007, 10**-1.55, r'$\chi^2$ = {0:0.3f}'.format(redchisqrM22), fontsize=30)
        plt.text(0.007, 10**-1.75, r'$p$ = {0:0.2e}'.format(p_val), fontsize=30)
        legend = ax.legend(loc='lower left', prop={'size':30}, labelspacing=0.35)
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
        #plt.show()
        #plt.savefig('C:/Users/Brendan/Desktop/test_format2/171_%ix_%iy_E.pdf' % (m2[m],l2[m]), format='pdf')
        #plt.savefig('C:/Users/Brendan/Desktop/171_slice2_double_optimize/171A_%ii_%ij.jpeg' % (l,m))
        #plt.savefig('C:/Users/Brendan/Desktop/171_points_square/pixel_%ii_%ij_new.jpeg' % (l2[m],m2[m]))
        #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
        #plt.close()