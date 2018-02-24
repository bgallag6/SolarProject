# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 17:57:54 2018

@author: Brendan
"""

import numpy as np
import sys

# set variables from command line
directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])
size = int(sys.argv[4])

#directory = 'S:'
#date = '20130626'
#wavelength = 171
#size = 16


cube_temp = []

# load derotated cube chunks
for i in range(size):
    temp = np.load('%s/DATA/Temp/%s/%i/chunk_%i_of_%i.npy' % (directory, date, wavelength, i, size))
    cube_temp.append(temp)
    
cube_final = np.vstack(cube_temp)  # stack chunks into final derotated array

np.save('%s/DATA/Temp/%s/%i/derotated.npy'% (directory, date, wavelength), cube_final)