# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 16:12:55 2017

@author: Brendan
"""

import numpy as np



# load original array 
original = np.load('C:/Users/Brendan/Desktop/testing_171/spectra_1seg_smoothed.npy')
print original.shape
orig_shape = np.array([original.shape[0], original.shape[1], original.shape[2]])

# create memory-mapped array with similar datatype and shape to original array
mmap_arr = np.memmap('C:/Users/Brendan/Desktop/testing_171/spectra_mmap_1seg_smoothed.npy', dtype='float64', mode='w+', shape=(original.shape[0],original.shape[1],original.shape[2]))

# write data to memory-mapped array
mmap_arr[:] = original[:]

# flush memory changes to disk, then remove memory-mapped object
del mmap_arr

np.save('C:/Users/Brendan/Desktop/testing_171/spectra_mmap_shape_1seg_smoothed.npy', orig_shape)