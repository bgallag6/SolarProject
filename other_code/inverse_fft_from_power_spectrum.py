# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 19:03:59 2017

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



num_freq = 299
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # 12 second cadence for the others
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

a0 = 10**-10.
n0 = 2.5
c0 = 10**-5.
p0 = 0.
l0 = 0.
w0 = 0.

spectra = a0*freqs**(-n0)+c0+p0*np.exp(-0.5*(((np.log(freqs))-(l0))/w0)**2)

plt.figure()
plt.loglog(freqs, spectra)

spectra2 = np.zeros((599))
spectra2[1:300] = spectra
spectra2[300:599] = np.flipud(spectra)

#plt.loglog(freqs,spectra)

v = np.fft.ifft(spectra2)
v2 = np.real(v)

plt.figure()
plt.plot(v)

## determine frequency values that FFT will evaluate
time_step = 12  # add as argument in function call, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]

sig = v
sig_fft = fftpack.fft(sig)
powers = np.abs(sig_fft)[pidxs]
norm = len(sig)  # to normalize the power
#powers = ((powers/norm)**2)*(1./(sig.std()**2))*2  # without normalization gives same result !!

plt.figure()
plt.loglog(freqs, powers)