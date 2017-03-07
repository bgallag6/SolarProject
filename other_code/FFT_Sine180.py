# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 07:13:58 2016

@author: Brendan
"""
##########################################
### Check 180-second Sine Wave Results ###
##########################################

import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
#from gatspy.periodic import LombScargleFast
#from gatspy.periodic import LombScargle
from scipy import fftpack


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


y = y3+y5+y6


### Generate sine wave WITHOUT errors
signal_base = np.sin(omega * t_interp)
signalA = np.sin(omega * t_interp + phi)
signalB = np.sin(omega*2 * t_interp + phi)
signal1 = signalA + signalB

### Generate sine wave WITH errors
signalA_e = np.sin(omega * t_interp + phi) + errors
signalB_e = np.sin(omega*2 * t_interp + phi) + errors
signal2 = signalA_e + signalB_e

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

"""
### Gatspy, using the same frequencies as FFT
fmin = freqss[0]
fmax = freqss[len(freqss)-1]
N = len(freqss)
df = (fmax - fmin) / (N-1)
model4 = LombScargle().fit(t_interp, sig)
power4 = model4.score_frequency_grid(fmin,df,N)
model5 = LombScargleFast().fit(t_interp, sig)
power5 = model5.score_frequency_grid(fmin,df,N)
freqs2 = fmin + df * np.arange(N)

#np.savetxt('C:/Users/Brendan/Desktop/PHYS 326/FFT_powers_errors.txt', powers)
#np.savetxt('C:/Users/Brendan/Desktop/PHYS 326/LS_powers_errors.txt', power5)
#np.savetxt('C:/Users/Brendan/Desktop/PHYS 326/freqs.txt', freqs2)
"""

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
plt.savefig('C:/Users/Brendan/Desktop/sine_random.jpeg')
 
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
plt.ylim(0,0.3)
plt.vlines(1./90.,0,1.5, linestyles='dashed', linewidth=1.5)
plt.text(0.0113, 0.24, '0.0111 Hz = 90 s',fontsize=20)
plt.vlines(1./180.,0,1.5, linestyles='dashed', linewidth=1.5)
plt.text(0.006, 0.17, '0.00556 Hz = 180 s',fontsize=20)
plt.vlines(1./360.,0,1.5, linestyles='dashed', linewidth=1.5)
plt.text(0.0001, 0.12, '0.00278 Hz = 360 s',fontsize=20)

#plt.text(0.0058,1.05,'0.00556 Hz = 180 s',fontsize=20)
#plt.legend(loc='upper left', prop={'size':17})
plt.show()   
plt.savefig('C:/Users/Brendan/Desktop/sine_fft_random.jpeg')

  