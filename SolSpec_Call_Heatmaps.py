# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 00:00:48 2017

@author: Brendan
"""

### Example full-module function calls 

import numpy as np
import SolSpec as ss
#import h5py

"""
## generate heatmaps
"""

directory = 'F:/Users/Brendan/Desktop/SolarProject'
#directory = 'C:/Users/Brendan/Desktop/FileTree'
date = '20130815'
wavelength = 171

"""
import sys

directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])
"""
#HEATMAPS = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_param.npy')
#VISUAL = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_visual.npy')

ss.heatmap(directory= '%s' % (directory), date='%s' % (date), wavelength= wavelength)

#r = ss.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20130626', wavelength=193, path_name='C:/Users/Brendan/Desktop/test_delete/')