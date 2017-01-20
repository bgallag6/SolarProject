# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 18:44:53 2017

@author: Brendan
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

with h5py.File("F:/Users/Brendan/Desktop/SolarProject/hdf5/20141025_304A_2300_2600i_1700_2400j_float_dogbox.hdf5",'r') as f:

    h_map = np.array(f['params'])
        
    vis = np.array(f['visual'])

    R = np.zeros((h_map.shape[1],h_map.shape[2]))
    
    """
    #20130530 1600
    pla_min = -1.
    pla_max = 1.0e-7
    index_min = 1.75
    index_max = 2.5
    plc_min = -0.000004
    plc_max = 1.
    amp_min = 0.
    amp_max = 0.012
    loc_min = -10.
    loc_max = -5.75
    wid_min = 0.4
    wid_max = 0.8
    """
    
    #20141025 304
    pla_min = -1.
    pla_max = 4.0e-7
    index_min = 1.5
    index_max = 2.5
    plc_min = 0.001
    plc_max = 1.
    amp_min = 0.
    amp_max = 0.005
    loc_min = -6.5
    loc_max = -4.5
    wid_min = 0.
    wid_max = 0.55



for i in range(R.shape[0]):
    for j in range(R.shape[1]):
        R[i][j] = vis[0][i][j]
for i in range(R.shape[0]):
    for j in range(R.shape[1]):
        if h_map[0][i][j] > pla_max or h_map[0][i][j] < pla_min \
        or h_map[1][i][j] > index_max or h_map[1][i][j] < index_min \
        or h_map[2][i][j] > plc_max or h_map[2][i][j] < plc_min \
        or h_map[3][i][j] > amp_max or h_map[3][i][j] < amp_min \
        or h_map[4][i][j] > loc_max or h_map[4][i][j] < loc_min \
        or h_map[5][i][j] > wid_max or h_map[5][i][j] < wid_min:
            R[i][j] = np.nan
                
plt.imshow(R)