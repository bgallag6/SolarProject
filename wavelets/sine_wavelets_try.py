# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 14:47:56 2017

@author: Brendan
"""

import numpy as np
from waveletFunctions import wavelet, wave_signif
import matplotlib.pylab as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
#from gatspy.periodic import LombScargleFast
#from gatspy.periodic import LombScargle
from scipy import fftpack

__author__ = 'Evgeniya Predybaylo'


# WAVETEST Example Python script for WAVELET, using NINO3 SST dataset
#
# See "http://paos.colorado.edu/research/wavelets/"
# The Matlab code written January 1998 by C. Torrence is modified to Python by Evgeniya Predybaylo, December 2014
#
# Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
#   changed all "log" to "log2", changed logarithmic axis on GWS to
#   a normal axis.
# ------------------------------------------------------------------------------------------------------------------

# READ THE DATA
#sst = np.loadtxt('sst_nino3.dat')  # input SST time series

### Specify properties for sine wave
period = 180 # seconds
frequency = 1.0 / period # Hertz
omega = 2. * np.pi * frequency # radians per second
phi = 0.5 * np.pi # radians

### Create time steps at which to evaluate sine wave
t_interp = [(6*k) for k in range(0,360)]
t_interp = np.array(t_interp)
t_interp = t_interp.astype(float)

l = len(t_interp)

### Generate errors to add to sine wave
errors = scipy.stats.norm.rvs(loc=0, scale=.1, size=360)

ran = np.random.random((l))
ran2 = np.random.random((l))
ran3 = np.random.random((l))
ran4 = np.random.random((l))


y3 = 1.*np.sin((t_interp*omega/2.)+(np.pi/2.))+2.*ran
y5 = 1.5*np.sin((t_interp*omega))+1.5*ran2
y6 = 2.*np.sin((t_interp*omega*2.)-(np.pi/2.))+1.*ran3


### Generate sine wave WITHOUT errors
signal_base = np.sin(omega * t_interp)
signalA = np.sin(omega * t_interp + phi)
signalB = np.sin(omega*2 * t_interp + phi)
signal1 = signalA + signalB

### Generate sine wave WITH errors
signalA_e = np.sin(omega * t_interp + phi) + errors
signalB_e = np.sin(omega*2 * t_interp + phi) + errors
signal2 = signalA_e + signalB_e

#y = y3+y5+y6
y = signal_base

#plt.plot(t_interp, signal1)

#np.savetxt('C:/Users/Brendan/Desktop/PHYS 326/sin_pure.txt', signal1)
#np.savetxt('C:/Users/Brendan/Desktop/PHYS 326/sin_noise.txt', signal2)
           
### FFT
sig = y

time_step = 6
sample_freq = fftpack.fftfreq(sig.size, d=time_step)
sig_fft = fftpack.fft(sig)
pidxs = np.where(sample_freq > 0)
freqss = sample_freq[pidxs]
powers = np.abs(sig_fft)[pidxs]
norm = len(t_interp)
powers *= 2 / (norm * sig.std() ** 2)  # to normalize the power
p = fftpack.fft(sig)

#np.savetxt('C:/Users/Brendan/Desktop/PHYS 326/sin_pure.txt', signal1)
#np.savetxt('C:/Users/Brendan/Desktop/PHYS 326/sin_noise.txt', signal2)
#np.savetxt('C:/Users/Brendan/Desktop/PHYS 326/frequencies.txt', sample_freq)


### PLot   
fig = plt.figure(figsize=(20,11))
#plt.xlim(100,500)
#plt.xlim(10**-3, 10**-1.5)
#plt.ylim(0,1)
#plt.plot(1.0/freqss, powers/2, label='FFT')  # we divided the FFT power by 2 (seemed to normalize pretty close to gatspy)
#plt.plot(t_interp, sig, label='signal')  # we divided the FFT power by 2 (seemed to normalize pretty close to gatspy)
plt.plot(t_interp, y, linewidth=2.0, label='signal')  # we divided the FFT power by 2 (seemed to normalize pretty close to gatspy)
#plt.plot(1.0/freqs2, power4, 'ko', label='Gatspy')
#plt.plot(1.0/freqs2, power5, label='Gatspy Fast')
plt.grid()
#plt.title('Sine Wave w/ %i-second Period' % period, fontsize=20, fontweight='bold', y = 1.01)
plt.title('Sine Wave w/ ???-second Periods', fontsize=20, fontweight='bold', y = 1.01)
plt.ylabel('Amplitude', fontsize=20)
plt.xlabel('Time [s]', fontsize=20)
ran = np.abs(np.max(y)-np.min(y))
plt.ylim(np.min(y)-(0.2*ran), np.max(y)+(0.2*ran))
##plt.legend(loc='upper left', prop={'size':17})
plt.show()  
#plt.savefig('C:/Users/Brendan/Desktop/sine_random.jpeg')
 
### PLot   
fig = plt.figure(figsize=(20,11))
#plt.xlim(100,500)
#plt.xlim(10**-3, 10**-1.5)
#plt.ylim(0,1)
#plt.plot(1.0/freqss, powers/2, label='FFT')  # we divided the FFT power by 2 (seemed to normalize pretty close to gatspy)
#plt.plot(freqss[4:7], powers[4:7]/2., 'r', label='FFT')  # we divided the FFT power by 2 (seemed to normalize pretty close to gatspy)
#plt.plot(freqss[10:13], powers[10:13]/2., 'g', label='FFT')  # we divided the FFT power by 2 (seemed to normalize pretty close to gatspy)
plt.plot(freqss, powers/2., linewidth=2.0, label='FFT')  # we divided the FFT power by 2 (seemed to normalize pretty close to gatspy)
#plt.plot(1.0/freqs2, power4, 'ko', label='Gatspy')
#plt.plot(1.0/freqs2, power5, label='Gatspy Fast')
plt.grid()
plt.title('Fast Fourier Transform of Sine Wave w/ ???-second Periods', fontsize=20, fontweight='bold', y =1.01)
plt.ylabel('Power', fontsize=20, labelpad = 10)
plt.xlabel('Frequency [Hz]', fontsize=20)
plt.xlim(0,0.014)
plt.ylim(0,1.1)
plt.vlines(1./90.,0,1.5, linestyles='dashed', linewidth=1.5)
plt.text(0.0113, 0.24, '0.0111 Hz = 90 s',fontsize=20)
plt.vlines(1./180.,0,1.5, linestyles='dashed', linewidth=1.5)
plt.text(0.006, 0.17, '0.00556 Hz = 180 s',fontsize=20)
plt.vlines(1./360.,0,1.5, linestyles='dashed', linewidth=1.5)
plt.text(0.0001, 0.12, '0.00278 Hz = 360 s',fontsize=20)

#plt.text(0.0058,1.05,'0.00556 Hz = 180 s',fontsize=20)
#plt.legend(loc='upper left', prop={'size':17})
plt.show()   
#plt.savefig('C:/Users/Brendan/Desktop/sine_fft_random.jpeg')



sst = y

#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"
variance = np.std(sst, ddof=1) ** 2
sst = (sst - np.mean(sst)) / np.std(sst, ddof=1)
n = len(sst)
dt = 6/60.
#time = np.arange(len(sst)) * dt + 1871.0  # construct time array
time = t_interp/60.
#xlim = ([1870, 2000])  # plotting range
xlim = ([0,2154/60.])
pad = 1  # pad the time series with zeroes (recommended)
dj = 0.25  # this will do 4 sub-octaves per octave
s0 = 2 * dt  # this says start at a scale of 6 months
j1 = 7 / dj  # this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.72  # lag-1 autocorrelation for red noise background
#mother = 'PAUL'
mother = 'MORLET'

# Wavelet transform:
wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)
power = (np.abs(wave)) ** 2  # compute wavelet power spectrum

# Significance levels: (variance=1 for the normalized SST)
signif = wave_signif(([1.0]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother)
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand signif --> (J+1)x(N) array
sig95 = power / sig95  # where ratio > 1, power is significant

# Global wavelet spectrum & significance levels:
global_ws = variance * (np.sum(power, axis=1) / n)  # time-average over all times
dof = n - scale  # the -scale corrects for padding at edges
global_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=1, lag1=lag1, dof=dof, mother=mother)

# Scale-average between El Nino periods of 2--8 years
avg = np.logical_and(scale >= 1., scale < 10.)
Cdelta = 0.776  # this is for the MORLET wavelet
scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
scale_avg = power / scale_avg  # [Eqn(24)]
scale_avg = variance * dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2, lag1=lag1, dof=([2, 7.9]), mother=mother)

#------------------------------------------------------ Plotting

#--- Plot time series
plt.figure(figsize=(18, 9))
plt.subplot(221)
plt.plot(time, sst)
plt.xlim(xlim[:])
plt.xlabel('Time (minutes)')
plt.ylabel('Amplitude')
plt.title(r'a) Sine Wave | y = sin(2$\pi$/180)')
plt.hold(False)

# --- Plot 2--8 yr scale-average time series
plt.subplot(222)
plt.plot(time, scale_avg)
plt.xlim(xlim[:])
plt.xlabel('Time (minutes)')
plt.ylabel('Avg variance (amplitude^2)?')
plt.title('d) 1-10 min Scale-average Time Series')
plt.hold(True)
plt.plot(xlim, scaleavg_signif + [0, 0], '--')
plt.hold(False)

#--- Contour plot wavelet power spectrum
plt3 = plt.subplot(223)
levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
CS = plt.contourf(time, period, np.log2(power), len(levels))  #*** or use 'contour'
im = plt.contourf(CS, levels=np.log2(levels))
plt.xlabel('Time (minutes)')
plt.ylabel('Period (minutes)')
plt.title('b) Sine Wave Wavelet Power Spectrum (in base 2 logarithm)')
plt.xlim(xlim[:])
# 95# significance contour, levels at -99 (fake) and 1 (95# signif)
plt.hold(True)
plt.contour(time, period, sig95, [-99, 1], colors='k')
# cone-of-influence, anything "below" is dubious
plt.plot(time, coi, 'k')
plt.hold(False)
# format y-scale
plt3.set_yscale('log', basey=2, subsy=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt3.ticklabel_format(axis='y', style='plain')
plt3.invert_yaxis()
# set up the size and location of the colorbar
divider = make_axes_locatable(plt3)
cax = divider.append_axes("bottom", size="5%", pad=0.5)
plt.colorbar(im, cax=cax, orientation='horizontal')

#--- Plot global wavelet spectrum
plt4 = plt.subplot(224)
plt.plot(global_ws, period)
plt.hold(True)
plt.plot(global_signif, period, '--')
plt.hold(False)
plt.xlabel('Power (amplitude^2)')
plt.ylabel('Period (minutes)')
plt.title('c) Global Wavelet Spectrum')
plt.xlim([0, 1.25 * np.max(global_ws)])
# format y-scale
plt4.set_yscale('log', basey=2, subsy=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt4.ticklabel_format(axis='y', style='plain')
plt4.invert_yaxis()

plt.tight_layout()

plt.show()

#plt.savefig('C:/Users/Brendan/Desktop/sine_180sec_wavelet.pdf', format='pdf', bbox_inches='tight')

# end of code