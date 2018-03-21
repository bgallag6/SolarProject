# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 16:18:06 2017

@author: Brendan
"""

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import interp1d
import matplotlib


#DATA = np.load('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/derotated.npy')
#DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/data_sort/20130626/20130626_193_-450_-200i_-200_200j_data_rebin1.npy')
DATA = np.load('C:/Users/Brendan/Desktop/Files/solar/time_series_171.npy')
print(DATA.shape[0])

#Ex =  np.load('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/exposure.npy')
Ex =  np.load('F:/DATA/Temp/20130626/171/exposure.npy')   
#Ex = np.load('C:/Users/Brendan/Desktop/exposure.npy')

#TIME = np.load('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/time.npy')
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/data_sort/20130626/20130626_193_-450_-200i_-200_200j_time.npy') 
TIME = np.load('C:/Users/Brendan/Desktop/Files/solar/time_arr_171.npy')

#timeseries = np.load('C:/Users/Brendan/Desktop/timeseries_arr.npy')     
    
ex_interp = np.interp(TIME,Ex,Ex)  # interpolate pixel-intensity values onto specified time grid
  

m2 = [188, 726, 722, 872]
l2 = [523, 328, 1427, 875] # = 1, 3, 6, 24 of full DATA array
match = [1,3,6,24]

m2_title = ['(c) Power Law Dominated w/o Lorentzian', '(d) Power Law Dominated w/ Lorentzian', '(a) Tail Dominated w/o Lorentzian', '(b) Tail Dominated w/o Lorentzian']

point_label = ['C', 'D', 'A', 'B']

pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values

matplotlib.rc('text', usetex = True)  # use with latex commands
plt.rc('font',**{'family':'serif','serif':['Times']})

for l in range(1):
    
    for m in range(len(m2)): 
        
        
        #v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
        
        plt.rcParams["font.family"] = "Times New Roman"
        font_size = 27        
        
        # Plot models + display combined-model parameters + uncertainties
        #"""
        fig = plt.figure(figsize=(12,10))
        ax = plt.gca()
        plt.title('%s: Point %s' % (m2_title[m], point_label[m]), y = 1.01, fontsize=font_size)
        #plt.plot(t_interp,v_interp,'k')
        ax.tick_params(axis='both', which='major', pad=10)
        plt.xlim(0,720)
        plt.plot(TIME/60., DATA[match[m]]/ex_interp, 'k')
        plt.xlabel(r'Time [min]', fontsize=font_size, labelpad=10)
        plt.ylabel('Normalized Intensity', fontsize=font_size, labelpad=10)
        plt.xticks([0,120,240,360,480,600,720], fontsize=font_size)
        plt.yticks(fontsize=font_size)
        #plt.text(600, 8.5, r'$n$ = {0:0.2f}'.format(2.05), fontsize=font_size)
        #plt.text(600, 8.0, r'$\beta$ = {0:0.1f} [min]'.format(2.04), fontsize=font_size)
        #plt.text(600, 7.5, r'$\delta$ = {0:0.3f}'.format(2.03), fontsize=font_size)
        #plt.text(600, 7.0, r'$p$ = {0:0.3g}'.format(2.02), fontsize=font_size)
        #plt.text(600, 6.5, r'$r$ = {0:0.3g}'.format(2.01), fontsize=font_size)
        #plt.text(600, 7.5, r'$\alpha$ = {0:0.2f}'.format(1.04), fontsize=font_size)  # works here but not spectra points?

        #plt.savefig('C:/Users/Brendan/Desktop/sample_timeseries_%s_lorentz.pdf' % (point_label[m]), format='pdf', bbox_inches='tight')
        #plt.close()
        #"""