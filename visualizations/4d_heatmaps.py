# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 16:27:19 2017

@author: Brendan
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20101208'
wavelength = 1600

H = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
#H = np.load('%s/DATA/Output/%s/PCB/%i/param.npy' % (directory, date, wavelength))

titles = [r'Power Law Slope-Coefficient', r'Power Law Index', r'Power Law Tail', r'Gaussian Amplitude', r'Gaussian Location [min]', r'Gaussian Width', 'F-Statistic', r'Gaussian Amplitude Scaled', 'p-Value']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']

p1 = 5
p2 = 2

param1 = H[p1,50:-50,70:-70] # for 20160520
param2 = H[p2,50:-50,70:-70] # for 20160520

#param1 = H[p1,180:-175,350:-80] # for 20140818
#param2 = H[p2,180:-175,350:-80] # for 20140818
#param1 = H[p1,70:-55, 190:270] # for 20160426
#param2 = H[p2,70:-55, 190:270] # for 20160426
#param1 = H[p1,160:-130, 260:-165] # for 20160520
#param2 = H[p2,160:-130, 260:-165] # for 20160520
#param1 = H[p1,180:400,390:625] # for 20131118
#param2 = H[p2,180:400,390:625] # for 20131118
#param1 = H[p1,105:205,110:210] # for 20140606
#param2 = H[p2,105:205,110:210] # for 20140606
#param1 = H[p1,60:150,60:150] # for 20140112
#param2 = H[p2,60:150,60:150] # for 20140112

#param1 = H[p1,150:600,150:650] # for 20170329
#param2 = H[p2,150:600,150:650] # for 20170329
#param1 = H[p1,375:475,625:750] # for 20141025
#param2 = H[p2,375:475,625:750] # for 20141025
#param1 = H[p1,170:320,170:320] # for 20130626
#param2 = H[p2,170:320,170:320] # for 20130626


if p1 == 4:
    param1 = (1./np.exp(param1)) / 60.
if p2 == 4:
    param2 = (1./np.exp(param2)) / 60.
    cmap = 'jet_r'
else:
    cmap = 'jet'

fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111, projection='3d')

X = np.linspace(0,param1.shape[1],param1.shape[1])
Y = np.linspace(0,param2.shape[0],param2.shape[0])
X,Y = np.meshgrid(X,Y)

cbar = ax.scatter(X, Y, param1, c=param2, cmap=cmap)  #use vmin/vmax when necessary
plt.title('%s %i$\AA$: %s [z] vs %s [color]' % (date, wavelength, titles[p1], titles[p2]), y=1.1, fontsize=23)
ax.set_xlabel('x-position (pixels)')
ax.set_ylabel('y-position (pixels)')
ax.set_zlabel('%s' % titles[p1])
#ax.set_zlim(0,0.04)
fig.colorbar(cbar)

plt.show()