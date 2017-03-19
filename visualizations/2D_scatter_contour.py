# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 00:15:19 2017

@author: Brendan
"""

from plotly.graph_objs import *
import numpy as np

HEATMAPS = np.load('C:/Users/Brendan/Desktop/SDO/param_20130530_1600_2300_2600i_2200_3000j_data_rebin4b.npy')
hmap = HEATMAPS[0]

x = np.reshape(h_map[3], (h_map.shape[1]*h_map.shape[2]))
y = np.reshape(h_map[4], (h_map.shape[1]*h_map.shape[2]))

plot([Histogram2dContour(x=x, y=y, contours=Contours(coloring='heatmap')),
       Scatter(x=x, y=y, mode='markers', marker=Marker(color='white', size=3, opacity=0.3))], show_link=False)