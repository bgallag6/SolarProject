# -*- coding: utf-8 -*-
"""
Created on Wed May 03 16:41:19 2017

@author: Brendan
"""

import numpy as np
import plotly.plotly as py
from plotly.graph_objs import Surface
import plotly
import plotly.offline as offline
import plotly.graph_objs as go
from plotly.graph_objs import *
from plotly import tools
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
plotly.offline.init_notebook_mode() # run at the start of every ipython notebook
"""
directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20140818'
wavelength = 1600

H = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
#H = np.load('%s/DATA/Output/%s/PCB/%i/param.npy' % (directory, date, wavelength))

titles = [r'Power Law Slope-Coefficient', r'Power Law Index', r'Power Law Tail', r'Gaussian Amplitude', r'Gaussian Location [min]', r'Gaussian Width', 'F-Statistic', r'Gaussian Amplitude Scaled', 'p-Value']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']

p1 = 5
p2 = 2

param1 = H[p1,180:-175,350:-80] # for 20140818
param2 = H[p2,180:-175,350:-80] # for 20140818

param1 = [[1,1,1],
          [2,2,2],
          [3,3,3]]


if p1 == 4:
    param1 = (1./np.exp(param1)) / 60.
if p2 == 4:
    param2 = (1./np.exp(param2)) / 60.
    cmap = 'jet_r'
else:
    cmap = 'jet'

z1 = [param1
]
z2 = z1
z3 = z1


plot([
    dict(z=z1, type='surface'),
    dict(z=z2, showscale=False, opacity=0.9, type='surface'),
    dict(z=z3, showscale=False, opacity=0.9, type='surface')])
"""    

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20121018'
wavelength1 = 171

H_1600 = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength1))
#H = np.load('%s/DATA/Output/%s/PCB/%i/param.npy' % (directory, date, wavelength))

titles = [r'Power Law Slope-Coefficient', r'Power Law Index', r'Power Law Tail', r'Gaussian Amplitude', r'Gaussian Location [min]', r'Gaussian Width', 'F-Statistic', r'Gaussian Amplitude Scaled', 'p-Value']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']

p1 = 1
p2 = 2

H_1600 = H_1600[:]


#param1 = H_1600[p1,0:,0:] # for 20160520
#param2 = H[p2,0:,0:] # for 20160520

#param1 = H[p1,180:-175,350:-80] # for 20140818
#param2 = H[p2,180:-175,350:-80] # for 20140818
#param1 = H[p1,30:160,40:170] # for 20140819
#param2 = H[p2,30:160,40:170] # for 20140819
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
#param1 = H[p1, 325:-250, 350:-125] # for 20160414
#param2 = H[p2, 325:-250, 350:-125] # for 20160414

#param1 = H[p1,150:600,150:650] # for 20170329
#param2 = H[p2,150:600,150:650] # for 20170329
#param1 = H[p1,375:475,625:750] # for 20141025
#param2 = H[p2,375:475,625:750] # for 20141025
#param1 = H[p1,170:320,170:320] # for 20130626
#param2 = H[p2,170:320,170:320] # for 20130626

h1600 = np.ndarray.tolist(H_1600)


p0_1600 = Surface(
    z=h1600[0],
    visible=False,
)

p1_1600 = Surface(
    z=h1600[1],
) 

p2_1600 = Surface(
    z=h1600[2],
    visible=False,
)

p3_1600 = Surface(
    z=h1600[3],
    visible=False,
)

p4_1600 = Surface(
    z=h1600[4],
    visible=False,
)

p5_1600 = Surface(
    z=h1600[5],
    visible=False,
)

p6_1600 = Surface(
    z=h1600[6],
    visible=False,
)

p7_1600 = Surface(
    z=h1600[7],
    visible=False,
)

    
data = Data([p0_1600, p1_1600, p2_1600, p3_1600, p4_1600, p5_1600, p6_1600, p7_1600])
layout = Layout(
    title='Parameter 4D Surface Plots',
    annotations=[dict(text='1600A',
                      font=dict(size=18, color='Black'),
                      x=-0.15, y=1.0,
                      xref='paper', yref='paper',
                      showarrow=False)],
    updatemenus=list([
        dict(
            x=-0.05,
            y=0.95,
            yanchor='top',
            buttons=list([
                dict(
                    args=['visible', [True, False, False, False, False, False, False, False]],
                    label='Slope Coefficient',
                    method='restyle',
                ),
                dict(
                    args=['visible', [False, True, False, False, False, False, False, False]],
                    label='Power Law Index',
                    method='restyle',
                ),
                dict(
                    args=['visible', [False, False, True, False, False, False, False, False]],
                    label='Power Law Tail',
                    method='restyle',
                ),
                dict(
                    args=['visible', [False, False, False, True, False, False, False, False]],
                    label='Gaussian Amplitude',
                    method='restyle',
                ),
                dict(
                    args=['visible', [False, False, False, False, True, False, False, False]],
                    label='Gaussian Location',
                    method='restyle',
                ),
                dict(
                    args=['visible', [False, False, False, False, False, True, False, False]],
                    label='Gaussian Width',
                    method='restyle',
                ),
                dict(
                    args=['visible', [False, False, False, False, False, False, True, False]],
                    label='F-Statistic',
                    method='restyle',
                ),
                dict(
                    args=['visible', [False, False, False, False, False, False, False, True]],
                    label='Gaussian Amplitude Scaled',
                    method='restyle',
                )
            ]),
        ),
    ]),
)
fig = Figure(data=data, layout=layout)
plot(fig)