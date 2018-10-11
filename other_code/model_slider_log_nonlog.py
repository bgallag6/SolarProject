# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 10:49:54 2018

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy import fftpack 

 
def LorentzianM2(f, A, n, C, P, fp, fw):
    return A*f**-n + C + P*(1./ (1.+((f-fp)/fw)**2))

def Lorentzian(f, P, fp, fw):
    return P*(1./ (1.+((f-fp)/fw)**2))

def LorentzianM2B(f, A, n, C, P, fp, fw):
    return A*f**-n + C + P*(1./ (1.+((np.log(f)-fp)/fw)**2))

def LorentzianB(f, P, fp, fw):
    return P*(1./ (1.+((np.log(f)-fp)/fw)**2))

    
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
fp0 = (1./300.)
fp0B = -5.5
fw0 = 0.001
fw0B = 0.1
#s = a0*np.sin(2*np.pi*f0*t)
if model == "lorentzian":
    s = LorentzianM2(f,A0,n0,C0,P0,fp0,fw0)
    s0 = Lorentzian(f,P0,fp0,fw0)
    sB = LorentzianM2B(f,A0,n0,C0,P0,fp0B,fw0B)
    s0B = LorentzianB(f,P0,fp0B,fw0B)


# setup plot with initial parameters
fig, ax = plt.subplots(figsize=(15,10))
plt.subplots_adjust(bottom=0.4)
#l, = plt.plot(f, s, lw=2, color='red')
#l0, = plt.plot(f, s0, lw=2, color='green')
#lB, = plt.plot(f, sB, lw=2, color='red', linestyle='dashed')
#l0B, = plt.plot(f, s0B, lw=2, color='green', linestyle='dashed')
l, = plt.loglog(f, s, lw=2, color='red', label=r'$\nu$')
l0, = plt.loglog(f, s0, lw=2, color='green', label=r'$\ln(\nu)$')
lB, = plt.loglog(f, sB, lw=2, color='red', linestyle='dashed')
l0B, = plt.loglog(f, s0B, lw=2, color='green', linestyle='dashed')
plt.xlim(10**-4.5, 10**-1)
plt.ylim(10**-6, 10**0)
plt.vlines((1./180.),10**-6,10**0, linestyle='dotted')
plt.vlines((1./300.),10**-6,10**0, linestyle='dashed')
plt.vlines((1./600.),10**-6,10**0, linestyle='dashed')
plt.vlines((1./1200.),10**-6,10**0, linestyle='dashed')
plt.legend(fontsize=17)

# setup parameter value sliders
axcoeff = plt.axes([0.2, 0.26, 0.6, 0.02])
axindex = plt.axes([0.2, 0.22, 0.6, 0.02])
axamp = plt.axes([0.2, 0.18, 0.6, 0.02])
axloc = plt.axes([0.2, 0.14, 0.6, 0.02])
axlocB = plt.axes([0.2, 0.10, 0.6, 0.02])
axwid = plt.axes([0.2, 0.06, 0.6, 0.02])
axwidB = plt.axes([0.2, 0.02, 0.6, 0.02])

scoeff = Slider(axcoeff, 'Coeff', 1e-10, 0.0001, valinit=A0)
sindex = Slider(axindex, 'Index', 0.3, 4.0, valinit=n0)
samp = Slider(axamp, 'Amp', 0.00001, 1., valinit=P0)
sloc = Slider(axloc, 'Loc', (1./660.), (1./60.), valinit=fp0)
slocB = Slider(axlocB, 'LocB', -6.5, -4.6, valinit=fp0B)
swid = Slider(axwid, 'Wid', 0.0001, 0.1, valinit=fw0)
swidB = Slider(axwidB, 'WidB', 0.05, 0.8, valinit=fw0B)

resetax = plt.axes([0.8, 0.02, 0.1, 0.02])
button = Button(resetax, 'Reset', hovercolor='0.975')


def update(val):
    A2 = scoeff.val
    n2 = sindex.val
    C2 = C0
    P2 = samp.val
    fp2 = sloc.val
    fp2B = slocB.val
    fw2 = swid.val
    fw2B = swidB.val
    #l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    print('%0.3e' % P2, '%0.3f' % fp2, '%0.3f' % fw2)
    if model == "lorentzian":
        s2 = LorentzianM2(f,A2,n2,C2,P2,fp2,fw2)
        s20 = Lorentzian(f,P2,fp2,fw2)
        s2B = LorentzianM2B(f,A2,n2,C2,P2,fp2B,fw2B)
        s20B = LorentzianB(f,P2,fp2B,fw2B)
    l.set_ydata(s2)
    l0.set_ydata(s20)
    lB.set_ydata(s2B)
    l0B.set_ydata(s20B)
    fig.canvas.draw_idle()
scoeff.on_changed(update)
sindex.on_changed(update)
samp.on_changed(update)
sloc.on_changed(update)
slocB.on_changed(update)
swid.on_changed(update)
swidB.on_changed(update)


def reset(event):
    scoeff.reset()
    sindex.reset()
    samp.reset()
    sloc.reset()
    slocB.reset()
    swid.reset()
    swidB.reset()
button.on_clicked(reset)

plt.show()