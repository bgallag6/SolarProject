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

#DATA = np.load('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/derotated.npy')
#DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/data_sort/20130626/20130626_193_-450_-200i_-200_200j_data_rebin1.npy')
DATA = np.load('C:/Users/Brendan/Desktop/20160905/193_0_cube.npy')
print DATA.shape[0]

#Ex =  np.load('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/exposure.npy')
#Ex =  np.load('F:/Users/Brendan/Desktop/SolarProject/data_sort/20130626/20130626_193_-450_-200i_-200_200j_exposure.npy')   
Ex = np.load('C:/Users/Brendan/Desktop/20160905/exposure.npy')

#TIME = np.load('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/time.npy')
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/data_sort/20130626/20130626_193_-450_-200i_-200_200j_time.npy') 
TIME = np.load('C:/Users/Brendan/Desktop/20160905/time.npy')

#timeseries = np.load('C:/Users/Brendan/Desktop/timeseries_arr.npy')     
    
t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters
  
#x2 = [722, 867, 765, 757, 525, 708, 743]
#y2 = [1441, 864, 325, 319, 551, 352, 322]

#x2 = [1]
#y2 = [1]

#x2 = [722, 525, 757, 743]
#y2 = [1441, 551, 319, 322]

m2_title = ['Tail Dominated w/o Gaussian', 'Power-Law Dominated w/o Gaussian', 'Power-Law Dominated w/ Gaussian', 'Tail Dominated w/ Gaussian']

pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values

for l in range(1):
    
    #for m in range(len(x2)): 
    for m in range(1):
            
        for k in range(0,DATA.shape[0]):
              #im=DATA[k]/(Ex[k])	  # get image + normalize by exposure time  (time went nuts?)
              im=DATA[k]
              #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])  # finds pixel-box median
              #pixmed[k]= im[y2[m],x2[m]]	# median  <-- use this
              pixmed[k]= im[m]	# median  <-- use this

        pixmed = pixmed/Ex  # normalize by exposure time    
        
        # The derotation introduces some bad data towards the end of the sequence. This trims that off
        bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
        last_good_pos = bad - 1			# retain only data before the <=zero
        
        # Get time and pixel values
        #v=pixmed[0:last_good_pos]		
        #t=TIME[0:last_good_pos]
        v=pixmed		
        t=TIME
    
        v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
 
        # Plot models + display combined-model parameters + uncertainties
        #"""
        fig = plt.figure(figsize=(15,15))
        #plt.title('%s: Pixel %ix, %iy' % (m2_title[m], x2[m], y2[m]), y = 1.01, fontsize=30)
        #plt.title('Pixel (%i, %i)' % (x2[m], y2[m]), y = 1.01, fontsize=25)
        plt.plot(t_interp,v_interp,'k')
        plt.xlabel('Time [Seconds]', fontsize=25, labelpad=10)
        plt.ylabel('Intensity', fontsize=25, labelpad=10)
        plt.xticks(fontsize=19)
        plt.yticks(fontsize=19)
        #plt.show()
        #plt.savefig('/mnt/data-solar/Gallagher/171_timeseries_%ix_%iy.pdf' % (x2[m], y2[m]), format='pdf')
        #plt.savefig('C:/Users/Brendan/Desktop/test_format/171_%ix_%iy_timeseries.pdf' % (x2[m],y2[m]), format='pdf')
        #plt.savefig('C:/Users/Brendan/Desktop/171_%ix_%iy_timeseries.pdf' % (x2[m],y2[m]), format='pdf')
        #plt.savefig('C:/Users/Brendan/Desktop/title25_labels21_tick17.pdf', format='pdf')
        #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
        #plt.close()
        #"""