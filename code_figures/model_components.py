# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 11:44:15 2017

@author: Brendan
"""

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt   

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Power-Law-fitting function (Model M1)
def PowerLaw1(f, A, n):
    return A*f**-n
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)


freqs2 = np.linspace((1./7188.), (1./24.04), 299)
freqs3 = np.linspace((1./7188.), (1./24.04), 299)



fit1 = PowerLaw1(freqs2, 10**-8., 2.)  # could get rid of this if not making smaller m2_fit
fit = GaussPowerBase(freqs2, 10**-8., 2., 10**-4.2, 0.01, -5.5, 0.3)  # could get rid of this if not making smaller m2_fit
m2G_fit = Gauss(freqs2, 0.01, -5.5, 0.3)
m2P_fit = PowerLaw(freqs2, 10**-8., 2., 10**-4.2)

noise = np.random.normal(1,0.2,299)
noise_M = fit*noise


fig = plt.figure(figsize=(20,15))
ax = plt.gca()
#plt.title(r'Model Components: $A$ $x^{-n}$ + $C$ + $P$ $e^{{-\frac{(\lnx-\beta)^{2}}{\sigma^{2}}}}$', y = 1.01, fontsize=30)
plt.title(r'Model Components', y = 1.01, fontsize=30)
#$(C/A)^{-\frac{1}{n}}$
plt.ylim((10**-5.,10**0))
plt.xlim((10**-4.,10**-1.3))
plt.xticks(fontsize=21)
plt.yticks(fontsize=21)
ax.tick_params(axis='both', which='major', pad=5)

#plt.loglog(f,s,'k')
#plt.loglog(f_fit, m1_fit, label='Power Law - M1')
#plt.loglog(freqs2, m2P_fit, 'g', label='Power Law - M2')
#plt.loglog(freqs2, m2G_fit, 'g--', label='M2 - Gaussian', linewidth=1.5)
plt.loglog(freqs2, m2P_fit, 'g', label='Power Law', linewidth=2.5)
plt.loglog(freqs2, m2G_fit, 'g--', label='Gaussian', linewidth=2.5)
#plt.loglog(f_fit, m2_fit, 'r', label='Combined - M2')
plt.loglog(freqs2, fit, 'purple', label='Combined Model', linewidth=2.5)
#plt.loglog(freqs2, fit1, linewidth=1.5)
plt.loglog(freqs2, noise_M, 'k', linewidth=1.)

plt.xlabel('Frequency (Hz)', fontsize=21, labelpad=5)
plt.ylabel('Power', fontsize=21, labelpad=5)
#plt.hlines((fit1[0]), 10**-4., 10**-3.5, linestyles='dashed', label='Slope Coefficient', linewidth=2.)
#plt.loglog(freqs2[0:15],fit[0:15], label='Power-Law Index', linewidth=2.)
#plt.loglog(freqs2[75:],fit[75:], label='PowerLaw Tail', linewidth=2.)
#plt.vlines(10**-1.5, 10**-5., 10**-4.00, linestyles='dashed', color='grey', label='Power-Law Tail', linewidth=2.)
#plt.vlines((1.0/(np.exp(5.47))),10**-3.23,10**-1.99, linestyles='dashed', color='purple', label='Gaussian Amplitude', linewidth=2.)
#plt.vlines((1.0/(np.exp(5.5))),10**-8,10**-0.5, linestyles='dashed', color='red', label='Gaussian Location', linewidth=2.)
#plt.hlines(10**-4.5, 10**-2.83, 10**-1.93, linestyles='dashed', color='orange', label='Gaussian Width', linewidth=2.)

legend = ax.legend(loc='upper right', prop={'size':27}, labelspacing=0.35)
for label in legend.get_lines():
    label.set_linewidth(3.0)  # the legend line width
plt.show()
plt.savefig('C:/Users/Brendan/Desktop/model_components.jpeg')
#plt.savefig('C:/Users/Brendan/Desktop/model_components.pdf', format='pdf')
#plt.close()
