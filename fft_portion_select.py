# -*- coding: utf-8 -*-
"""
Created on Sun Jan 08 16:54:19 2017

@author: Brendan
"""

"""
# widgets example code: slider_demo.py
The SpanSelector is a mouse widget to select a xmin/xmax range and plot the
detail view of the selected region in the lower axes
*** make so that can view full timeseries, can select portion and see if 
*** something pops out from FFT
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

fig = plt.figure(figsize=(16, 12))
#ax = fig.add_subplot(211)
ax1 = plt.subplot2grid((3,5),(0, 0), colspan=5, rowspan=1)

#f = np.loadtxt('C:/Users/Brendan/Desktop/PHYS 326/older/seconds.txt')
#s = np.loadtxt('C:/Users/Brendan/Desktop/PHYS 326/older/values.txt')

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20120909'
wavelength = 171

derotated = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))
#num_freq = derotated.shape[2]/2  # determine nubmer of frequencies that are used
#freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script

t = np.array([12*i for i in range(3600)]) 

s = derotated[0][0]

ax1.plot(t, s, '-')
ax1.set_title('Timeseries -- (Select segment to compute FFT)', fontsize=20, y=1.01)
ax1.set_xlim(t.min(), t.max())
ax1.set_ylim(s.min(), s.max())

#ax2 = fig.add_subplot(212)
ax2 = plt.subplot2grid((3,5),(1, 0), colspan=5, rowspan=2)

if wavelength == 1600 or wavelength == 1700:
    #time_step = 24  # 24 second cadence for these wavelengths
    time_step = 12  # Jacks dataset
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
    m1_fit = PowerLaw(f, A, n, C)  
    line3.set_data(f,m1_fit)
    
   
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
