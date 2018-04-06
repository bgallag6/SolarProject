# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 20:57:55 2018

@author: Brendan
"""

"""
############################
############################
# Create & Save Memory-Mapped Array
############################
############################
"""

import numpy as np
import sys
import os


# set variables from command line
file_path = sys.argv[1]
delete_temp = sys.argv[2]

# load original array 
original = np.load('%s.npy' % file_path)


# determine dimensions of array
if original.ndim == 1:
    orig_shape = np.array([original.shape[0]])
    
elif original.ndim == 2:
    orig_shape = np.array([original.shape[0], original.shape[1]])
    
elif original.ndim == 3:
    orig_shape = np.array([original.shape[0], original.shape[1], original.shape[2]])
    
elif original.ndim == 4:
    orig_shape = np.array([original.shape[0], original.shape[1], original.shape[2], original.shape[3]])
    
   
# create memory-mapped array with similar datatype and shape to original array
mmap_arr = np.memmap('%s_mmap.npy' % file_path, dtype='%s' % original.dtype, mode='w+', shape=tuple(orig_shape))


# write data to memory-mapped array
mmap_arr[:] = original[:]


# flush memory changes to disk, then remove memory-mapped object and original array
del mmap_arr
del original


# save memory-mapped array dimensions to use when loading
np.save('%s_mmap_shape.npy' % file_path, orig_shape)


# delete original array if specified
if delete_temp == "y":
    if os.path.isfile('%s.npy' % file_path):
        os.remove('%s.npy' % file_path)
    else:    ## Show an error ##
        print("Error: %s.npy file not found" % file_path)