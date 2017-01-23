# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 23:36:57 2017

@author: Brendan
"""

import plotly.plotly as py
from plotly.graph_objs import *
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go

HEATMAPS = np.load('C:/Users/Brendan/Desktop/SDO/param_20130530_1600_2300_2600i_2200_3000j_data_rebin4b.npy')

h0 = HEATMAPS[0]
h1 = HEATMAPS[1]
h2 = HEATMAPS[2]
h3 = HEATMAPS[3]
h4 = HEATMAPS[4]
h5 = HEATMAPS[5]

trace1 = go.Contour(
    z=h0,
    line=dict(smoothing=0.85),
    visible=False,
    colorbar=dict(
            title='Slope Coefficient',
            titleside='right',
            titlefont=dict(
                size=14,
                family='Arial, sans-serif'))
)

trace2 = go.Contour(
    z=h1,
    line=dict(smoothing=0.85),
    colorbar=dict(
            title='Index Value',
            titleside='right',
            titlefont=dict(
                size=14,
                family='Arial, sans-serif'))
)

trace3 = go.Contour(
    z=h2,
    line=dict(smoothing=0.85),
    visible=False,
    colorbar=dict(
            title='Tail Value',
            titleside='right',
            titlefont=dict(
                size=14,
                family='Arial, sans-serif'))
)

trace4 = go.Contour(
    z=h3,
    line=dict(smoothing=0.85),
    visible=False,
    colorbar=dict(
            title='Gaussian Amplitude',
            titleside='right',
            titlefont=dict(
                size=14,
                family='Arial, sans-serif'))
)

trace5 = go.Contour(
    z=h4,
    line=dict(smoothing=0.85),
    visible=False,
    colorbar=dict(
            title='Gaussian Location',
            titleside='right',
            titlefont=dict(
                size=14,
                family='Arial, sans-serif'))   
)

trace6 = go.Contour(
    z=h5,
    line=dict(smoothing=0.85),
    visible=False,
    colorbar=dict(
            title='Gaussian Width',
            titleside='right',
            titlefont=dict(
                size=14,
                family='Arial, sans-serif'))    
)

data = Data([trace1, trace2, trace3, trace4, trace5, trace6])
layout = Layout(
    title='Parameter Heatmaps',
    updatemenus=list([
        dict(
            x=-0.05,
            y=1,
            yanchor='top',
            buttons=list([
                dict(
                    args=['visible', [True, False, False, False, False, False]],
                    label='Power Law Coefficient',
                    method='restyle'
                ),
                dict(
                    args=['visible', [False, True, False, False, False, False]],
                    label='Power Law Index',
                    method='restyle'
                ),
                dict(
                    args=['visible', [False, False, True, False, False, False]],
                    label='Power Law Tail',
                    method='restyle'
                ),
                dict(
                    args=['visible', [False, False, False, True, False, False]],
                    label='Gaussian Amplitude',
                    method='restyle'
                ),
                dict(
                    args=['visible', [False, False, False, False, True, False]],
                    label='Gaussian Location',
                    method='restyle'
                ),
                dict(
                    args=['visible', [False, False, False, False, False, True]],
                    label='Gaussian Width',
                    method='restyle'
                )
            ]),
        )
    ]),
)
fig = Figure(data=data, layout=layout)
plot(fig)