# -*- coding: utf-8 -*-
"""
Created on Sat Mar 04 09:10:50 2017

@author: Brendan
"""

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *

timeseries = np.load('C:/Users/Brendan/Desktop/timeseries_arr.npy')     
    
x2 = [722, 525, 757, 743]
y2 = [1441, 551, 319, 322]

point_label = ['A', 'B', 'C', 'D']

m2_title = ['Tail Dominated w/o Gaussian', 'Power-Law Dominated w/o Gaussian', 'Power-Law Dominated w/ Gaussian', 'Tail Dominated w/ Gaussian']

TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/data_sort/20130626/20130626_171_-500_500i_-500_500j_time.npy') 

timeseries = np.load('C:/Users/Brendan/Desktop/timeseries_arr.npy')     
    
t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters
t_interp = t_interp/60
  
for l in range(1):
    
    for m in range(len(x2)):    
            
        # Plot models + display combined-model parameters + uncertainties
        #"""
        fig = plt.figure(figsize=(15,15))
        #ax = plt.gca()  # get current axis -- to set colorbar 
        plt.title('%s: Point %s' % (m2_title[m], point_label[m]), y = 1.01, fontsize=30)
        #plt.title('Pixel (%i, %i)' % (x2[m], y2[m]), y = 1.01, fontsize=25)
        plt.plot(t_interp,timeseries[m],'k')
        plt.xlabel('Time [min]', fontsize=30, labelpad=10)
        plt.ylabel('Intensity', fontsize=30, labelpad=10)
        plt.xlim(0,np.max(t_interp))
        plt.xticks([0,120,240,360,480,600,720], fontsize=30)
        plt.yticks(fontsize=30)
        #ax.tick_params(axis='both', which='major', pad=15)
        #plt.show()
        #plt.savefig('/mnt/data-solar/Gallagher/171_timeseries_%ix_%iy.pdf' % (x2[m], y2[m]), format='pdf')
        plt.savefig('C:/Users/Brendan/Desktop/171_point_%s_timeseries.pdf' % (point_label[m]), format='pdf')
        #plt.savefig('C:/Users/Brendan/Desktop/171_%ix_%iy_timeseries.pdf' % (x2[m],y2[m]), format='pdf')
        #plt.savefig('C:/Users/Brendan/Desktop/title25_labels21_tick17.pdf', format='pdf')
        #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
        #plt.close()
        #"""
        