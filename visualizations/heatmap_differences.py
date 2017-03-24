# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 21:06:00 2017

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


p_6_double = np.load('C:/Users/Brendan/Desktop/20130626_171_test/171/param_6_double.npy')
p_6_single = np.load('C:/Users/Brendan/Desktop/20130626_171_test/171/param_6_single.npy')
p_4_double = np.load('C:/Users/Brendan/Desktop/20130626_171_test/171/param_4_double.npy')
p_4_single = np.load('C:/Users/Brendan/Desktop/20130626_171_test/171/param_4_single.npy')
#arth_double = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_param_slope6_arthm_dogboxTRF.npy')

#diff_double = np.abs(p_6_double - p_4_double) / p_6_double
diff_6 = p_6_double - p_6_single
diff_double = p_6_double - p_4_double
diff_6_double = p_6_double - p_4_single

arr = [diff_6, diff_double, diff_6_double]

titles = ['6/double vs 6/single', '6/double vs 4/double', '6/double vs 4/single']
titles_p = [r'Slope Coeff', r'Power Law Index', r'Power Law Tail', r'Gaussian Amplitude', r'Gaussian Location', r'Gaussian Width']
save_t = ['6_single','4_double','4_single']
save_title = ['coeff','index','tail','amp','loc','wid']
"""
for i in range(2):
    flat_param2 = np.reshape(diff_6[i], (diff_6[i].shape[0]*diff_6[i].shape[1]))
    flat_param3 = np.reshape(diff_double[i], (diff_6[i].shape[0]*diff_6[i].shape[1]))
    flat_param4 = np.reshape(diff_6_double[i], (diff_6[i].shape[0]*diff_6[i].shape[1]))

    fig = plt.figure(figsize=(12,9))
    plt.title('Spectra Fitting Methods Comparison: %s' % titles[i], y = 1.01, fontsize=25)
    #plt.xlabel('%s' % cbar_labels[i], fontsize=20, labelpad=10)
    plt.ylabel('Bin Count', fontsize=20, labelpad=10)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    #plt.xlim(M2_low[i],M2_high[i])
    y, x, _ = plt.hist(flat_param2, bins=200, color='black', label = '6 - Single')
    y, x, _ = plt.hist(flat_param3, bins=200, color='blue', alpha=0.5, label = '4 - Double')
    y, x, _ = plt.hist(flat_param4, bins=200, color='red', alpha=0.5, label = '4 - Single')
    plt.ylim(0, y.max()*1.1)
    plt.legend(loc = 'upper right')
    #plt.hist(flatten_slopes, bins='auto')  # try this (actually think we want constant bins throughout wavelengths)
    #plt.savefig('%s/%s_%i_Histogram_%s.jpeg' % (path_name, date, wavelength, names[i]))
    #plt.savefig('%s/%s_%i_Histogram_%s.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
"""
#"""
#diff_geo_arth = np.abs((geo - arth) / geo) * 100.
#diff_arth_arth2 = np.abs((arth - arth_double) / arth) * 100.
#diff_geo_arth2 = np.abs((geo - arth_double) / geo) * 100.

for i in range(6):
    for c in range(len(arr)):
        plt.figure(figsize=(12,10))
        plt.title('%s -- %s' % (titles[c], titles_p[i]), fontsize = 25, y=1.01)
        h_min = np.percentile(arr[c][i],5)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(arr[c][i],95)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        #plt.imshow(diff_geo_arth[i], vmin=np.min(diff_geo_arth[i]), vmax=np.max(diff_geo_arth[i]))
        #plt.imshow(np.flipud(diff_double[i]), vmin=0., vmax=100.)
        plt.imshow(np.flipud(arr[c][i]), vmin=h_min, vmax = h_max)
        plt.colorbar()
        #plt.savefig('C:/Users/Brendan/Desktop/171_diff/%s_%s.pdf' % (save_title[i],save_t[c]), format='pdf')
#"""