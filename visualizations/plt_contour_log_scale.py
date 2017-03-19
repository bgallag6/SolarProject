# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 22:53:34 2017

@author: Brendan
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal

HEATMAPS = np.load('C:/Users/Brendan/Desktop/SDO/param_20130530_1600_2300_2600i_2200_3000j_data_rebin4b.npy')
hmap = HEATMAPS[3]

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
z = hmap


# Automatic selection of levels works; setting the
# log locator tells contourf to use a log scale:
fig, ax = plt.subplots()
cs = ax.contourf(X, Y, z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)

# Alternatively, you can manually set the levels
# and the norm:
#lev_exp = np.arange(np.floor(np.log10(z.min())-1),
#                    np.ceil(np.log10(z.max())+1))
#levs = np.power(10, lev_exp)
#cs = P.contourf(X, Y, z, levs, norm=colors.LogNorm())

# The 'extend' kwarg does not work yet with a log scale.

cbar = fig.colorbar(cs)

plt.show()