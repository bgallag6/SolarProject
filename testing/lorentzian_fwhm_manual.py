# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 13:03:26 2018

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy import fftpack 
import scipy.signal

def GaussianM2(f, A, n, C, P, fp, fw):
    return A*f**-n + C + P*np.exp(-0.5*((np.log(f)-fp)/fw)**2)
    #return A*f**-n

#"""    
def LorentzianM2(f, A, n, C, P, fp, fw):
    #return A*f**-n + C + P*(1./ (1.+((np.log(f)-fp)/fw)**2))
    return A*f**-n + C + P*(1./ (1.+((np.log(f)-fp)/fw)**2))

def Lorentzian(f, P, fp, fw):
    #return A*f**-n + C + P*(1./ (1.+((np.log(f)-fp)/fw)**2))
    return P*(1./ (1.+((np.log(f)-fp)/fw)**2))
#"""

"""
def LorentzianM2(f, A, n, C, P, fp, fw, a):
    #return A*f**-n + C + P*(1./ (1.+((np.log(f)-fp)/fw)**2))
    return A*f**-n + C + P*(1./ (1.+((np.log(f)-fp)/(2*fw/(1+np.exp(a*(np.log(f)-fp)))))**2))
"""
    
#SPECTRA = np.load('C:/Users/Brendan/Desktop/SDO/spectra_20120923_211A_(528)_(132)x_(100)_100y.npy')
    
## load in array of segment-averaged pixel FFTs
#spectra_array = SPECTRA

#num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used

   
# determine frequencies that FFT will evaluate
num_freq = 99999
#num_freq = 299
freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
time_step = 12  # add as argument, or leave in as constant?
#time_step = 24  # add as argument, or leave in as constant?
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]
f = freqs

df = np.log10(freqs[1:len(freqs)]) - np.log10(freqs[0:len(freqs)-1])
df2 = np.zeros_like(freqs)
df2[0:len(df)] = df
df2[len(df2)-1] = df2[len(df2)-2]
ds = df2


model = "lorentzian"

arr = np.zeros((6,4))
arr2 = np.zeros((6,4))

for loc in range(2,8):
    for wid in [1,2,3,4]:
    #for wid in [4]:
        A0 = 0.000001
        n0 = 1.
        C0 = 0.0001
        P0 = 0.01
        #fp0 = -5.5
        #fp0 = np.log(1./(8.*60.))
        fp0 = np.log(1./(loc*60.))
        #fw0 = 0.1
        fw0 = wid*0.1
        a0 = 0.
        #s = a0*np.sin(2*np.pi*f0*t)
        if model == "lorentzian":
            s = LorentzianM2(f,A0,n0,C0,P0,fp0,fw0)
            #s = LorentzianM2(f,A0,n0,C0,P0,fp0,fw0,a0)
            #s += np.random.randn((num_freq))*s*0.1
        elif model == "gaussian":
            s = GaussianM2(f,A0,n0,C0,P0,fp0,fw0)
        
        
        """
        # setup plot with initial parameters
        fig, ax = plt.subplots(figsize=(8,8))
        #l, = plt.loglog(f, s, lw=2, color='red')
        l = plt.loglog(f, s, lw=2, color='red')
        #plt.xlim(10**-4.5, 10**-1)
        #plt.ylim(10**-6, 10**0)
        plt.vlines((1./180.),10**-6,10**0, linestyle='dotted')
        plt.vlines((1./300.),10**-6,10**0, linestyle='dashed')
        plt.vlines((1./600.),10**-6,10**0, linestyle='solid')
        plt.ylim(10**-6,10**1)
        """
        
        
        M2_low = [0., 0.3, -0.01, 0.00001, -6.5, 0.05]
        M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]                         
        #M2_low = [0., 0.3, -0.01, 0.00001, -6.5, 0.05, -100.]
        #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8, 100.]
        
        nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(LorentzianM2, f, s, p0=[A0,n0,C0,P0,fp0,fw0], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
        
        m2 = LorentzianM2(f, nlfit_gp[0],nlfit_gp[1],nlfit_gp[2],nlfit_gp[3],nlfit_gp[4],nlfit_gp[5])
        m2L = Lorentzian(f, nlfit_gp[3],nlfit_gp[4],nlfit_gp[5])
        #m2 = LorentzianM2(f, nlfit_gp[0],nlfit_gp[1],nlfit_gp[2],nlfit_gp[3],nlfit_gp[4],nlfit_gp[5],nlfit_gp[6])
        
        """
        #plt.loglog(f, m2)
        
        fig, ax = plt.subplots(figsize=(8,8))
        plt.vlines((1./180.),0.,1., linestyle='dotted')
        plt.vlines((1./300.),0.,1., linestyle='dashed')
        #plt.vlines((1./600.),0.,1., linestyle='solid')
        #plt.hlines(P0/2.,0.002,0.01)
        plt.plot(f, s, lw=2, color='red')
        plt.plot(f, m2L, 'g')
        plt.plot(f, m2)
        plt.ylim(0.,0.15)
        plt.xlim(0.,0.025)
        """
        
        """
        fig, ax = plt.subplots(figsize=(8,8))
        #plt.vlines((1./180.),0.,1., linestyle='dotted')
        #plt.vlines((1./300.),0.,1., linestyle='dashed')
        #plt.vlines((1./600.),0.,1., linestyle='solid')
        #plt.hlines(P0/2.,0.002,0.01)
        plt.title('Location = %i Minutes | Width = 0.%i' % (loc,wid))
        #plt.plot(f, s, lw=2, color='red')
        #plt.plot(f, m2L, 'g')
        #plt.plot(f, m2)
        #plt.ylim(0.,0.02)
        #plt.xlim(0.,0.015)
        plt.loglog(f, s, lw=2, color='red')
        plt.loglog(f, m2L, 'g')
        plt.loglog(f, m2)
        plt.ylim(10**-4,10**-1)
        plt.xlim(10**-4.5,10**-1)
        """
        
        widthp = nlfit_gp[5]
        
        b_freq = np.exp(nlfit_gp[4])
        idx = (np.abs(f-b_freq)).argmin()
        val0 = f[idx]
        b_time = (1./b_freq)/60.
        
        freq_l = f[:idx]
        freq_r = f[idx:]
        m2_l = m2L[:idx]
        m2_r = m2L[idx:]
        freq_l_min = np.argmin(m2_l)
        m2_l0 = m2_l[freq_l_min:]
        f_l0 = freq_l[freq_l_min:]
        
        #print(nlfit_gp[3])
        value = nlfit_gp[3]/2.
        idx_l = (np.abs(m2_l0-value)).argmin()
        val_l = f_l0[idx_l]
        
        idx_r = (np.abs(m2_r-value)).argmin()
        val_r = freq_r[idx_r]
        
        """
        plt.hlines(nlfit_gp[3]/2.,val_l,val_r)
        plt.vlines(val_l,0.,0.15)
        plt.vlines(val_r,0.,0.15)
        """
        
        print(widthp, (1./val_l) - (1./val_r), nlfit_gp[3])
        fwhm = (1./val_l) - (1./val_r)
        fwhm = 1. / (val_r - val_l)
        
        arr[loc-2][wid-1] = fwhm
        
        bleft = np.log(b_freq)-widthp
        bright = np.log(b_freq)+widthp
        #pleft = 1./(np.exp(bleft))
        #pright = 1./(np.exp(bright))
        pleft = (np.exp(bleft))
        pright = (np.exp(bright))
        
        #bleft = np.log10(b_freq)-(widthp/2)
        #bright = np.log10(b_freq)+(widthp/2)
        #pleft = 1./(10**bleft)
        #pright = 1./(10**bright)
        
        #arr2[loc-2][wid-1] = pleft-pright
        arr2[loc-2][wid-1] = 1. / (pright - pleft)
        """
        plt.figure()
        plt.plot([0.1,0.2,0.3,0.4],[97,211,334,444], 'r', label='8 min')
        plt.plot([0.1,0.2,0.3,0.4],[63,125,193,267], 'b', label='5 min')
        plt.plot([0.1,0.2,0.3,0.4],[45,95,147,201], 'g', label='3.75 min')
        plt.legend()
        
        points = [0.1,0.2,0.3,0.4]
        data = [[97,211,334,444],[63,125,193,267],[45,95,147,201]]
        
        def liner(x,m,b):
            return m*x + b
        
        for i in range(3):
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(liner, points, data[i])
            print(nlfit_gp[0])
        """
        #plt.savefig('C:/Users/Brendan/Desktop/loc_%i_wid_%i_loglog.pdf' % (loc,wid), format='pdf', bbox_inches='tight')