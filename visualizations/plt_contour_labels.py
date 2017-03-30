# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 22:40:34 2017

@author: Brendan
"""

"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

HEATMAPS = np.load('C:/Users/Brendan/Desktop/solar/results/20130626_final/20130626_193_-500_500i_-500_600j_param_slope6_arthm.npy')
hmap = HEATMAPS[1]

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

delta = 1.
x = np.arange(0., hmap.shape[1], delta)
y = np.arange(0., hmap.shape[0], delta)
X, Y = np.meshgrid(x, y)
#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
#Z = 10.0 * (Z2 - Z1)
Z = hmap

# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
plt.figure()
#CS = plt.contour(X, Y, Z)
CS = plt.contour(X, Y, Z, levels=[0.3,1.3])
plt.clabel(CS, inline=1, fontsize=10)
plt.colorbar()
plt.title('193A Contour : Power Law Index')