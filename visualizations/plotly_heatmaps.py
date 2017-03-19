# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 00:07:06 2017

@author: Brendan
"""

import plotly.plotly as py
from plotly.graph_objs import *
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go

HEATMAPS = np.load('C:/Users/Brendan/Desktop/SDO/param_20130530_193_2300_2600_2200_3000_float_numpy.npy')

h0 = HEATMAPS[0]
h1 = HEATMAPS[1]
h2 = HEATMAPS[2]
h3 = HEATMAPS[3]
h4 = HEATMAPS[4]
h5 = HEATMAPS[5]

trace1 = go.Heatmap(
    z=h0,
    visible=False,
)

trace2 = go.Heatmap(
    z=h1,
)

trace3 = go.Heatmap(
    z=h2,
    visible=False,
)

trace4 = go.Heatmap(
    z=h3,
    visible=False,
)

trace5 = go.Heatmap(
    z=h4,
    visible=False,
)

trace6 = go.Heatmap(
    z=h5,
    visible=False,
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