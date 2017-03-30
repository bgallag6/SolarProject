# -*- coding: utf-8 -*-
"""
Created on Sun Jan 08 16:54:19 2017

@author: Brendan
"""

"""
# widgets example code: slider_demo.py
The SpanSelector is a mouse widget to select a xmin/xmax range and plot the
detail view of the selected region in the lower axes
*** make so that can view full timeseries, can select portion and see if 
*** something pops out from FFT
"""

#"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from scipy import fftpack

fig = plt.figure(figsize=(16, 12))
#ax = fig.add_subplot(211)
ax1 = plt.subplot2grid((3,5),(0, 0), colspan=5, rowspan=1)

f = np.loadtxt('C:/Users/Brendan/Desktop/PHYS 326/older/seconds.txt')
s = np.loadtxt('C:/Users/Brendan/Desktop/PHYS 326/older/values.txt')

ax1.plot(f, s, '-')
ax1.set_title('Timeseries -- (Select segment to compute FFT)', fontsize=20, y=1.01)
ax1.set_xlim(f.min(), f.max())
ax1.set_ylim(s.min(), s.max())

#ax2 = fig.add_subplot(212)
ax2 = plt.subplot2grid((3,5),(1, 0), colspan=5, rowspan=2)

time_step = 12
sig = s
sample_freq = fftpack.fftfreq(sig.size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqss = sample_freq[pidxs]
sig_fft = fftpack.fft(sig)
powers = np.abs(sig_fft)[pidxs]
norm = len(sig)  # to normalize the power
powers = ((powers/norm)**2)*(1./(sig.std()**2))*2


ax2.set_xlim(10**-4.5, 10**-1)
ax2.set_ylim(10**-7, 10**0)        
line2, = ax2.loglog(freqss, powers, '-')


        
def onselect(fmin, fmax):
    indmin, indmax = np.searchsorted(f, (fmin, fmax))
    indmax = min(len(f) - 1, indmax)
    
    this_sig = s[indmin:indmax]
    this_sample_freq = fftpack.fftfreq(this_sig.size, d=time_step)
    this_pidxs = np.where(this_sample_freq > 0)
    this_freqss = this_sample_freq[this_pidxs]
    this_sig_fft = fftpack.fft(this_sig)
    this_powers = np.abs(this_sig_fft)[this_pidxs]
    this_norm = len(this_sig)  # to normalize the power
    this_powers = ((this_powers/this_norm)**2)*(1./(this_sig.std()**2))*2
    
    line2.set_data(this_freqss, this_powers)
    
    fig.canvas.draw()

# set useblit True on gtkagg for enhanced performance
span = SpanSelector(ax1, onselect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='red'))


plt.show()
#"""
