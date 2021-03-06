# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 16:04:57 2018

@author: Brendan
"""

import numpy as np
import sys
import yaml
import os

size = int(sys.argv[1])

stream = open('specFit_config.yaml', 'r')
cfg = yaml.load(stream)

directory = cfg['fits_dir']
date = cfg['date']
wavelength = cfg['wavelength']
mmap_derotate = cfg['mmap_derotate']
save_temp = cfg['save_temp']

#directory = 'S:'
#date = '20130626'
#wavelength = 171
#size = 16

cube_temp = []

# load derotated cube chunks
for i in range(size):
    temp = np.load('%s/DATA/Temp/%s/%i/chunk_%i_of_%i.npy' % (directory, date, wavelength, i+1, size))
    cube_temp.append(temp)
    
cube_final = np.vstack(cube_temp)  # stack chunks into final derotated array

del cube_temp


if mmap_derotate == "y":
    orig_shape = np.array([cube_final.shape[0], cube_final.shape[1], cube_final.shape[2]])
    
    # create memory-mapped array with similar datatype and shape to original array
    mmap_arr = np.memmap('%s/DATA/Temp/%s/%i/derotated_mmap.npy' % (directory, date, wavelength), dtype='%s' % cube_final.dtype, mode='w+', shape=tuple(orig_shape))
    
    # write data to memory-mapped array
    mmap_arr[:] = cube_final[:]
    
    # save memory-mapped array dimensions to use when loading
    np.save('%s/DATA/Temp/%s/%i/derotated_mmap_shape.npy' % (directory, date, wavelength), orig_shape)

    # save original array if specified
    if save_temp == "y":
        np.save('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength), cube_final)
        
    if save_temp == "n":
        for j in range(size):
            
            fn = '%s/DATA/Temp/%s/%i/chunk_%i_of_%i.npy' % (directory, date, wavelength, j+1, size)
            
            ## if file exists, delete it ##
            if os.path.isfile(fn):
                os.remove(fn)
            else:    ## Show an error ##
                print("Error: %s file not found" % fn)
    
    # flush memory changes to disk, then remove memory-mapped object and original array
    del mmap_arr
    del cube_final
    
else:
    np.save('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength), cube_final)