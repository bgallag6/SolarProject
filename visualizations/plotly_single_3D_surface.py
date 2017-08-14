# -*- coding: utf-8 -*-
"""
Created on Fri May 19 18:44:02 2017

@author: Brendan
"""
import matplotlib.pyplot as plt
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

directory = 'F:/Users/Brendan/Desktop/SolarProject'
#directory = 'S:'
date = '20140819'
wavelength1 = 1600

H_1600 = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength1))


titles = [r'Power Law Slope-Coefficient', r'Power Law Index', r'Power Law Tail', r'Gaussian Amplitude', r'Gaussian Location [min]', r'Gaussian Width', 'F-Statistic', r'Gaussian Amplitude Scaled', 'p-Value']
names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'p_value']

p1 = 1

#plt.imshow(H_304[1])

#H_1600 = H_1600[:,180:-175,350:-80]
#H_1600 = H_1600[:,100:-115,290:-20]

#param1 = H_1600[p1,0:,0:] # for 20160520
#param2 = H[p2,0:,0:] # for 20160520

#param1 = H[p1,180:-175,350:-80] # for 20140818
#param2 = H[p2,180:-175,350:-80] # for 20140818
#param1 = H[p1,30:160,40:170] # for 20140819
#param2 = H[p2,30:160,40:170] # for 20140819 -- 40:135,55:160
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

c = 5

if c == 4:
    H_1600[c] = (1./np.exp(H_1600[c])) / 60.

h1600 = np.ndarray.tolist(H_1600)




p0_1600 = Surface(
    z=h1600[0],
    #surfacecolor=h1600[c],
    #surfacecolor=h304[0],
    visible=False,
)

p1_1600 = Surface(
    z=h1600[1],
    #surfacecolor=h1600[c],
    #surfacecolor=h304[1]
) 

p2_1600 = Surface(
    z=h1600[2],
    #surfacecolor=h1600[c],
    #surfacecolor=h304[2],
    visible=False,
)

p3_1600 = Surface(
    z=h1600[3],
    #surfacecolor=h1600[c],
    #surfacecolor=h304[3],
    visible=False,
)

p4_1600 = Surface(
    z=h1600[4],
    #surfacecolor=h1600[c],
    #surfacecolor=h304[4],
    visible=False,
)

p5_1600 = Surface(
    z=h1600[5],
    #surfacecolor=h1600[c],
    #surfacecolor=h304[5],
    visible=False,
)

p6_1600 = Surface(
    z=h1600[6],
    #surfacecolor=h1600[c],
    #surfacecolor=h304[6],
    visible=False,
)

p7_1600 = Surface(
    z=h1600[7],
    #surfacecolor=h1600[c],
    #surfacecolor=h304[7],
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