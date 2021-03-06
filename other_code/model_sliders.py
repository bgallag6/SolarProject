# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 10:57:58 2017

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy import fftpack 

def GaussianM2(f, A, n, C, P, fp, fw):
    return A*f**-n + C + P*np.exp(-0.5*((np.log(f)-fp)/fw)**2)
    #return A*f**-n
    
def LorentzianM2(f, A, n, C, P, fp, fw):
    return A*f**-n + C + P*(1./ (1.+((np.log(f)-fp)/fw)**2))

    
#SPECTRA = np.load('C:/Users/Brendan/Desktop/SDO/spectra_20120923_211A_(528)_(132)x_(100)_100y.npy')
    
## load in array of segment-averaged pixel FFTs
#spectra_array = SPECTRA

#num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used

   
# determine frequencies that FFT will evaluate
num_freq = 299
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # add as argument, or leave in as constant?
#time_step = 24  # add as argument, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]
f = freqs


model = "lorentzian"

A0 = 0.0001
n0 = 1.
C0 = 0.0001
P0 = 0.1
fp0 = -5.5
fw0 = 0.1
#s = a0*np.sin(2*np.pi*f0*t)
if model == "lorentzian":
    s = LorentzianM2(f,A0,n0,C0,P0,fp0,fw0)
elif model == "gaussian":
    s = GaussianM2(f,A0,n0,C0,P0,fp0,fw0)


# setup plot with initial parameters
fig, ax = plt.subplots(figsize=(15,10))
plt.subplots_adjust(bottom=0.4)
l, = plt.loglog(f, s, lw=2, color='red')
plt.xlim(10**-4.5, 10**-1)
plt.ylim(10**-6, 10**0)
plt.vlines((1./180.),10**-6,10**0, linestyle='dotted')
plt.vlines((1./300.),10**-6,10**0, linestyle='dashed')
plt.vlines((1./600.),10**-6,10**0, linestyle='dashed')
plt.vlines((1./1200.),10**-6,10**0, linestyle='dashed')

# setup parameter value sliders
axcoeff = plt.axes([0.2, 0.26, 0.6, 0.02])
axindex = plt.axes([0.2, 0.22, 0.6, 0.02])
axtail = plt.axes([0.2, 0.18, 0.6, 0.02])
axamp = plt.axes([0.2, 0.14, 0.6, 0.02])
axloc = plt.axes([0.2, 0.1, 0.6, 0.02])
axwid = plt.axes([0.2, 0.06, 0.6, 0.02])

scoeff = Slider(axcoeff, 'Coeff', 1e-10, 0.0001, valinit=A0)
sindex = Slider(axindex, 'Index', 0.3, 4.0, valinit=n0)
stail = Slider(axtail, 'Tail', 0.000001, 0.03, valinit=C0)
samp = Slider(axamp, 'Amp', 0.00001, 1., valinit=P0)
sloc = Slider(axloc, 'Loc', -6.5, -4.6, valinit=fp0)
swid = Slider(axwid, 'Wid', 0.05, 0.8, valinit=fw0)

resetax = plt.axes([0.8, 0.02, 0.1, 0.02])
button = Button(resetax, 'Reset', hovercolor='0.975')


def update(val):
    A2 = scoeff.val
    n2 = sindex.val
    C2 = stail.val
    P2 = samp.val
    fp2 = sloc.val
    fw2 = swid.val
    #l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    if model == "lorentzian":
        s2 = LorentzianM2(f,A2,n2,C2,P2,fp2,fw2)
    elif model == "gaussian":
        s2 = GaussianM2(f,A2,n2,C2,P2,fp2,fw2)
    l.set_ydata(s2)
    fig.canvas.draw_idle()
scoeff.on_changed(update)
sindex.on_changed(update)
stail.on_changed(update)
samp.on_changed(update)
sloc.on_changed(update)
swid.on_changed(update)


def reset(event):
    scoeff.reset()
    sindex.reset()
    stail.reset()
    samp.reset()
    sloc.reset()
    swid.reset()
button.on_clicked(reset)

plt.show()