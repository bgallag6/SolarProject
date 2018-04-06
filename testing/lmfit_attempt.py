# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 13:55:03 2018

@author: Brendan
"""

from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, Model
import numpy as np
from scipy import fftpack  # doesnt work in module when called here???
import matplotlib.pyplot as plt

def GaussPowerBase(t, A2, n2, C2, P2, fp2, fw2):
    return A2*t**-n2 + C2 + P2*(1./ ((np.pi*fw2)*(1.+((np.log(t)-fp2)/fw2)**2)))

directory = 'F:'
date = '20130626'
wavelength = 171
    
cube_shape = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength))
spectra_array = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))

## load in array of segment-averaged pixel FFTs
SPECTRA = spectra_array
stddev = np.memmap('%s/DATA/Temp/%s/%i/uncertainties_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))

num_freq = 299

    
# determine frequency values that FFT will evaluate
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
if wavelength == 1600:
    time_step = 24
else:
    time_step = 12  # add as argument, or leave in as constant?

#time_step = 24  # add as argument, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]


m2 = [188, 726, 722, 872] # 727, 722, 189, 187, 867, 765, 757, 525, 708, 743, 790, 790, 790, 794, 796, 797, 798, 858, 861, 863, 872, 872, 876]
l2 = [523, 328, 1427, 875] # used in paper
#l2 = [524, 328, 1428, 875]


m2_title = ['(c) Power Law Dominated w/o Lorentzian', '(d) Power Law Dominated w/ Lorentzian', '(a) Tail Dominated w/o Lorentzian', '(b) Tail Dominated w/o Lorentzian']
#m2_title = 
point_label = ['C', 'D', 'A', 'B']


#for l in range(321,322):
for l in range(1):
    
    #for m in range(0,10):
    #for m in range(900,901):
    for m in range(len(m2)): 
    #for m in range(17,18):
        
                                        
        t = freqs  # frequencies
        s = spectra_array[l2[m]][m2[m]]  # fourier power
        #s = spectra_array[l2][m2]  # fourier power
        
        params = Parameters()
        params.add('A2',  value = 10**-7, min=0., max=0.01)
        params.add('n2',  value = 1.67, min=0.3, max=5.0)
        params.add('C2',  value = 0., min=-0.01, max=0.01)
        params.add('P2',  value = 0.001, min=10**-6, max=0.1)
        params.add('fp2', value = -5.5, min=-6.5, max=-4.6)
        params.add('fw2', value = 0.425, min=0.05, max=0.8)
        
        model = Model(GaussPowerBase, independent_vars=['t'])
        
        df = np.log10(t[1:len(t)]) - np.log10(t[0:len(t)-1])
        df2 = np.zeros_like(t)
        df2[0:len(df)] = df
        df2[len(df2)-1] = df2[len(df2)-2]
        ds = df2**-1.
        
        result = model.fit(s, params, t=t, weights=ds)
        
        report_fit(result)
        print(result.values['n2'])
        
        plt.figure()
        plt.loglog(t, s)  # data
        plt.loglog(t, GaussPowerBase(t=t, **result.values))  # best-fit model
        plt.loglog(t,ds)
        plt.xlim(10**-4, 10**-1)
        plt.ylim(10**-4, 10**0)