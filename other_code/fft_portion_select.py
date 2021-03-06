# -*- coding: utf-8 -*-
"""
Created on Sun Jan 08 16:54:19 2017

@author: Brendan
"""

"""
* make so that can view full timeseries, can select portion and see if something pops out from FFT
"""

#"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from scipy import fftpack

import numpy as np
import scipy.signal
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
from scipy import fftpack
from astropy.convolution import convolve, Box1DKernel
from numpy.random import randn
from mpi4py import MPI
from scipy.stats.stats import pearsonr   

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)

fig = plt.figure(figsize=(16, 12))
#ax = fig.add_subplot(211)
ax1 = plt.subplot2grid((3,5),(0, 0), colspan=5, rowspan=1)

#f = np.loadtxt('C:/Users/Brendan/Desktop/PHYS 326/older/seconds.txt')
#s = np.loadtxt('C:/Users/Brendan/Desktop/PHYS 326/older/values.txt')

#directory = 'F:/Users/Brendan/Desktop/SolarProject'
directory = 'S:'
date = '20140822'
wavelength = 1600

derotated = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))
#num_freq = derotated.shape[2]/2  # determine nubmer of frequencies that are used
#freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script

#t = np.array([12*i for i in range(3600)])
#t = np.array([24*i for i in range(1797)])  
t = np.array([24*i for i in range(derotated.shape[0])])  

#s = derotated[200][380]
#s = derotated[:,202,380]
s = derotated[:,50,100]

ax1.plot(t, s, '-')
ax1.set_title('Timeseries -- (Select segment to compute FFT)', fontsize=20, y=1.01)
ax1.set_xlim(t.min(), t.max())
ax1.set_ylim(s.min(), s.max())

#ax2 = fig.add_subplot(212)
ax2 = plt.subplot2grid((3,5),(1, 0), colspan=5, rowspan=2)

if wavelength == 1600 or wavelength == 1700:
    time_step = 24  # 24 second cadence for these wavelengths
    #time_step = 12  # Jacks dataset
else:
    time_step = 12  # 12 second cadence for the others
    #time_step = 24  # for half-cadence test
sig = s
sample_freq = fftpack.fftfreq(sig.size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]
sig_fft = fftpack.fft(sig)
powers = np.abs(sig_fft)[pidxs]
norm = len(sig)  # to normalize the power
powers = ((powers/norm)**2)*(1./(sig.std()**2))*2


ax2.set_xlim(10**-4.5, 10**-1)
ax2.set_ylim(10**-7, 10**0)      
ax2.vlines((1./180.),10**-7,10**0,linestyle='dotted')  
ax2.vlines((1./300.),10**-7,10**0,linestyle='dashed') 
line2, = ax2.loglog(freqs, powers, '-')
line3, = ax2.loglog(0)


        
def onselect(tmin, tmax):
    indmin, indmax = np.searchsorted(t, (tmin, tmax))
    indmax = min(len(t) - 1, indmax)

    
    this_sig = s[indmin:indmax]
    this_sample_freq = fftpack.fftfreq(this_sig.size, d=time_step)
    this_pidxs = np.where(this_sample_freq > 0)
    this_freqss = this_sample_freq[this_pidxs]
    this_sig_fft = fftpack.fft(this_sig)
    this_powers = np.abs(this_sig_fft)[this_pidxs]
    this_norm = len(this_sig)  # to normalize the power
    this_powers = ((this_powers/this_norm)**2)*(1./(this_sig.std()**2))*2
    
    line2.set_data(this_freqss, this_powers)
    
    f = this_freqss   
    
    # assign equal weights to all parts of the curve
    df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
    df2 = np.zeros_like(f)
    df2[0:len(df)] = df
    df2[len(df2)-1] = df2[len(df2)-2]
    ds = df2
    
    M1_low = [-0.002, 0.3, -0.01]
    M1_high = [0.002, 6., 0.01]
    nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, this_freqss, this_powers, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
    A, n, C = nlfit_l  # unpack fitting parameters
                                  
    M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
    M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]   

    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, this_freqss, this_powers, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays    
    
    A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
    
    nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, this_freqss, this_powers, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
               
    A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
    
    m2_param = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
                   
    # create model functions from fitted parameters
    #m1_fit = PowerLaw(f, A, n, C)        
    #m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
    m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22)
    line3.set_data(f,m2_fit2)
    
   
    #fig.canvas.draw()
    #A3 = ax2.text(0.0075, 10**-0.5, r'$A$ = {0:0.2e}'.format(A), fontsize=30)
    #n3 = ax2.text(0.0076, 10**-0.90, r'$n$ = {0:0.2f}'.format(n), fontsize=30)
    #plt.text(0.008, 10**-0.75, r'$C$ =  {0:0.3e}'.format(m2_param[2]), fontsize=25)
    #plt.text(0.007, 10**-0.73, r'$(C/A)^{-\frac{1}{n}}$ = %i [s]' % (1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1])), fontsize=30)
    #R3 = ax2.text(0.0075, 10**-1.35, r'$R$ = %1.0f [min]' % ((1./(C / A)**(-1./ n))/60.), fontsize=30)
    
    fig.canvas.draw()
    
    #ax2.texts.remove(A3)
    #ax2.texts.remove(n3)
    #ax2.texts.remove(R3)

# set useblit True on gtkagg for enhanced performance
span = SpanSelector(ax1, onselect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='red'))


plt.show()
#"""
