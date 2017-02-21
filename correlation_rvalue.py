# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 21:18:22 2017

@author: Brendan
"""
import numpy as np
import matplotlib.pyplot as plt

m2 = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/1600/20130626_1600_-450_-200i_-200_200j_m2.npy')
spectra = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/1600/20130626_1600_-450_-200i_-200_200j_spectra.npy')


from scipy.stats.stats import pearsonr

r = np.zeros((m2.shape[0], m2.shape[1]))

for i in range(m2.shape[0]):
    for j in range(m2.shape[1]):
        r_temp = pearsonr(m2[i][j], spectra[i][j])
        r[i][j] = r_temp[0]
        
plt.imshow(r)
