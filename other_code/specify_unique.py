# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 18:44:53 2017

@author: Brendan
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

#with h5py.File("F:/Users/Brendan/Desktop/SolarProject/hdf5/20130530_193A_2300_2600i_2200_3000j_float_dogbox.hdf5",'r') as f:

    #h_map = np.array(f['params'])
    #vis = np.array(f['visual'])

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20130626'
wavelength = 171

h_map = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
vis = np.load('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength))

#R = np.zeros((h_map.shape[1],h_map.shape[2]))

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
"""

"""
#20130530 193
pla_min = -1.
pla_max = 1.
index_min = 0.
index_max = 1.
plc_min = 0.0001
plc_max = 1.
amp_min = 0.
amp_max = 0.007
loc_min = -6.5
loc_max = -4.5
wid_min = 0.
wid_max = 0.8
"""

"""
#20130626 304 Filament
pla_min = -1.
pla_max = 2.0e-6
index_min = 1.3
index_max = 4.
plc_min = 0.0007
plc_max = 1.
amp_min = 0.
amp_max = 0.01
loc_min = -6.5
loc_max = -4.5
wid_min = 0.
wid_max = 0.8
"""

#"""
for k in range(1):
    #20130626 193 Coronal Hole
    pla_min = -1.
    pla_max = 1.
    index_min = 0.
    index_max = 1.5
    #index_max = 1.1 + (0.02*k)
    plc_min = 0.0018
    plc_max = 1.
    amp_min = 0.
    amp_max = 0.01
    loc_min = -6.5
    loc_max = -4.5
    wid_min = 0.
    wid_max = 0.8
    i_m = str(index_max)
#"""

    """
    #20130626 304 Coronal Hole
    for k in range(20):
        
        pla_min = -1.
        pla_max = 1.
        index_min = 0.
        #index_max = 1.3  # or 1.5, 1.4
        index_max = 1.0+(0.02*k)
        plc_min = 0.001
        plc_max = 1.
        amp_min = 0.
        amp_max = 0.01
        loc_min = -6.5
        loc_max = -4.5
        wid_min = 0.
        wid_max = 0.8
        
        i_m = str(index_max)
    """


    
    
    #for i in range(R.shape[0]):
    #    for j in range(R.shape[1]):
    #        R[i][j] = vis[0][i][j]
    R = vis[0,:-2,:-2]
    
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            if h_map[0][i][j] > pla_max or h_map[0][i][j] < pla_min \
            or h_map[1][i][j] > index_max or h_map[1][i][j] < index_min \
            or h_map[2][i][j] > plc_max or h_map[2][i][j] < plc_min \
            or h_map[3][i][j] > amp_max or h_map[3][i][j] < amp_min \
            or h_map[4][i][j] > loc_max or h_map[4][i][j] < loc_min \
            or h_map[5][i][j] > wid_max or h_map[5][i][j] < wid_min:
                R[i][j] = np.nan
                
    fig = plt.figure(figsize=(12,15))  
    #fig = plt.figure(figsize=(15,9))               
    plt.imshow(np.flipud(R))
    #plt.title('304A: Index Max = %s' % i_m)
    #plt.savefig('C:/Users/Brendan/Desktop/feature_detection/index_max_%s.jpeg' % i_m)

"""
## weighting attempt
R_binary = np.copy(R)

R_binary[R_binary > 0] = 1.

R_weight = np.copy(R_binary)

fig = plt.figure(figsize=(12,15))  
#fig = plt.figure(figsize=(15,9))               
plt.imshow(R_binary)

for i in range(1,R_binary.shape[0]-1):
    for j in range(1,R_binary.shape[1]-1):
        
        if R_binary[i][j] == 1:
            temp = np.zeros((8))
            temp[0] = R_binary[i-1][j-1]
            temp[1] = R_binary[i-1][j]
            temp[2] = R_binary[i-1][j+1]
            temp[3] = R_binary[i][j-1]
            temp[4] = R_binary[i][j+1]
            temp[5] = R_binary[i+1][j-1]
            temp[6] = R_binary[i+1][j]
            temp[7] = R_binary[i+1][j+1]
            
            R_weight[i][j] = np.sum(temp)

fig = plt.figure(figsize=(12,15))  
#fig = plt.figure(figsize=(15,9))               
plt.imshow(R_weight)      

R_weight2 = np.copy(R_weight)

R_binary2 = np.copy(R_weight)

R_binary2[R_binary2 > 0] = 1.

for i in range(1,R_binary2.shape[0]-1):
    for j in range(1,R_binary2.shape[1]-1):
        
        if R_binary2[i][j] == 1:
            temp = np.zeros((8))
            temp[0] = R_binary2[i-1][j-1]
            temp[1] = R_binary2[i-1][j]
            temp[2] = R_binary2[i-1][j+1]
            temp[3] = R_binary2[i][j-1]
            temp[4] = R_binary2[i][j+1]
            temp[5] = R_binary2[i+1][j-1]
            temp[6] = R_binary2[i+1][j]
            temp[7] = R_binary2[i+1][j+1]
            
            R_weight2[i][j] = np.sum(temp)

fig = plt.figure(figsize=(12,15))  
#fig = plt.figure(figsize=(15,9))               
plt.imshow(R_weight2)     
""" 