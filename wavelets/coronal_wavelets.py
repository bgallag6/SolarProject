# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 16:59:53 2017

@author: Brendan
"""

import numpy as np
from waveletFunctions import wavelet, wave_signif
import matplotlib.pylab as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import fftpack
import scipy.signal
import scipy.misc
from scipy import fftpack
from mpi4py import MPI
from scipy.stats.stats import pearsonr   
from scipy.stats import f as ff

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)

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
directory = 'S:'
date = '20140822'  # was 20130815, 171 for other figures
wavelength = 1600

data0 = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))
time0 = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))
exposure0 = np.load('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength))

img = data0[data0.shape[0]/2]

for x in range(90,100):

    #x,y = 105,142
    y = 141
    
    data_norm = data0[:,x,y] / exposure0
    sst = data_norm
    
    #----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------
    
    # normalize by standard deviation (not necessary, but makes it easier
    # to compare with plot on Interactive Wavelet page, at
    # "http://paos.colorado.edu/research/wavelets/plot/"
    variance = np.std(sst, ddof=1) ** 2
    sst = (sst - np.mean(sst)) / np.std(sst, ddof=1)
    n = len(sst)
    dt = 24/60.
    #time = np.arange(len(sst)) * dt + 1871.0  # construct time array
    time = time0/60.
    #xlim = ([1870, 2000])  # plotting range
    xlim = ([0,43200/60.])
    pad = 1  # pad the time series with zeroes (recommended)
    dj = 0.25  # this will do 4 sub-octaves per octave
    s0 = 2 * dt  # this says start at a scale of 6 months
    j1 = 7 / dj  # this says do 7 powers-of-two with dj sub-octaves each
    lag1 = 0.72  # lag-1 autocorrelation for red noise background
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
    avg = np.logical_and(scale >= 120/60., scale < 600/60.)
    Cdelta = 0.776  # this is for the MORLET wavelet
    scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
    scale_avg = power / scale_avg  # [Eqn(24)]
    scale_avg = variance * dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
    scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2, lag1=lag1, dof=([2, 7.9]), mother=mother)
    
    #------------------------------------------------------ Plotting
    
    #--- Plot time series
    plt.figure(figsize=(18, 11))
    plt.subplot(221)
    plt.plot(time, data_norm)
    plt.xlim(xlim[:])
    plt.xlabel('Time (minutes)')
    plt.ylabel('Intensity')
    plt.title('a) Intensity Time Series')
    plt.hold(False)
    
    ## determine frequency values that FFT will evaluate
    time_step = 24
      
    t_interp = np.linspace(0, time0[len(time0)-1], (time0[len(time0)-1]/time_step)+1)  # interpolate onto default-cadence time-grid
      
    n_segments = 6  # break data into 12 segments of equal length
    n = len(t_interp)
    rem = n % n_segments
    freq_size = (n - rem) / n_segments 
    
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    pixmed=np.empty(data0.shape[0])  # Initialize array to hold median pixel values
    spectra_seg = np.zeros((data0.shape[1],data0.shape[2],len(freqs)))
    
    
    pixmed = sst 
    
    v_interp = np.interp(t_interp,time0,pixmed)  # interpolate pixel-intensity values onto specified time grid
    
    data = v_interp
    
    avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
    
    data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
    split = np.split(data, n_segments)  # create split array for each segment
       
    for i in range(0,n_segments):               
        
      ## perform Fast Fourier Transform on each segment       
      sig = split[i]
      sig_fft = fftpack.fft(sig)
      #sig_fft = fftpack.rfft(sig)  # real-FFT                
      #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
      #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
      powers = np.abs(sig_fft)[pidxs]
      norm = len(sig)
      powers = ((powers/norm)**2)*(1./(sig.std()**2))*2   # normalize the power
      avg_array += powers
    
    avg_array /= n_segments
    
    f = freqs
    s = avg_array
    
    # assign equal weights to all parts of the curve
    df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
    df2 = np.zeros_like(f)
    df2[0:len(df)] = df
    df2[len(df2)-1] = df2[len(df2)-2]
    ds = df2
    #ds = 0.1*s
    
                                           
    ### fit data to models using SciPy's Levenberg-Marquart method
    
    """
    # initial guesses for fitting parameters
    M1_low = [-0.002, 0.3, -0.01]
    M1_high = [0.002, 6., 0.01]
    nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
    
      
    A, n, C = nlfit_l  # unpack fitting parameters
    
    # unpack uncertainties in fitting parameters from diagonal of covariance matrix
    dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
    """
    
    ## possibly use previous pixel's parameters as initial guesses for current pixel (issues creating wierd banding in images)
    
    ## fit data to combined power law plus gaussian component model
    #"""        
                                   
    M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
    M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
    #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
    
    # change method to 'dogbox' and increase max number of function evaluations to 3000
    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.43], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
    
    A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
    
    # unpack uncertainties in fitting parameters from diagonal of covariance matrix
    dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
    
    
      
    
    # create model functions from fitted parameters
    #m1_fit = PowerLaw(f_fit, A, n, C)
    m1_fit = PowerLaw(f, A2, n2, C2)
    amp_scale = PowerLaw(np.exp(fp2), A2, n2, C2)  # to extract the gaussian-amplitude scaling factor
    #m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
    m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
    #s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
    
    #"""
    nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
    #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
    A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
    dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
    m2_param = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
    m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22) 
    uncertainties = dA22, dn22, dC22, dP22, dfp22, dfw22  # do we want to keep a global array of uncertainties?
    
    residsM22 = (s - m2_fit2)
    chisqrM22 = ((residsM22/ds)**2).sum()
    redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6)  
    
    m2P_fit = PowerLaw(f, A22, n22, C22)  # only need if plotting
    m2G_fit = Gauss(f, P22, fp22, fw22)  # only need if plotting           
    #"""        
    
    #residsM2 = (s - s_fit_gp_full)
    residsM2 = (s - m2_fit)
    chisqrM2 = ((residsM2/ds)**2).sum()
    redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
    
    residsM1 = (s - m1_fit)
    chisqrM1 =  ((residsM1/ds)**2).sum()
    redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       
    
    # Plot models + display combined-model parameters + uncertainties
    residsM22 = (s - m2_fit2)
    chisqrM22 = ((residsM22/ds)**2).sum()
    redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6) 
        
    #f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))
    f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))
    
    #amp_scale = PowerLaw(np.exp(fp2), A2, n2, C2)  # to extract the gaussian-amplitude scaling factor
    amp_scale2 = PowerLaw(np.exp(fp22), A22, n22, C22)  # to extract the gaussian-amplitude scaling factor
    
    # generate p-value heatmap
    df1, df2 = 3, 6
    p_val = ff.sf(f_test2, df1, df2)
    
    r_val = pearsonr(m2_fit2, s)
    
    font_size = 17
    
    peak_fft = (1./np.exp(m2_param[4]))/60.
    
    plt.subplot(222)
    ax = plt.gca()
    plt.hold(True)
    plt.title('c) FFT Spectrum')
    plt.ylim((10**-4.5,10**-0.5))
    plt.xlim((10**-4.,10**-1.5))
    plt.loglog(f,s,'k', linewidth=1.5)
    plt.loglog(f, m2P_fit, 'g', label='M2 - Power Law', linewidth=1.5)
    plt.loglog(f, m2G_fit, 'g--', label='M2 - Gaussian', linewidth=1.5)
    plt.loglog(f, m2_fit2, 'purple', label='M2 - Combined', linewidth=1.5)
    #plt.loglog(f, m1_fit, linewidth=1.3)
    #plt.loglog(f, m2P_fit, 'g', label=r'$A\nu^{-n} + C$', linewidth=1.3)
    #plt.loglog(f, m2G_fit, 'g--', label=r'$\alpha\ e^{{-\frac{(\ln\nu-\beta)^{2}}{\sigma^{2}}}}$', linewidth=1.3)
    #plt.loglog(f, m2_fit2, 'purple', label='Combined Model', linewidth=1.3)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power')
    #plt.vlines((1.0/180.),10**-3.5,10**-0.5, linestyles='dotted', label='3 minutes')
    plt.vlines(np.exp(m2_param[4]),10**-4.5,10**-0.5, 'b', linestyles='dashed', label='%0.2f minutes' % peak_fft, linewidth=1.5)
    #plt.vlines((0.0093),10**-8,10**1, linestyles='dotted', label='3 minutes')
    
    #rect = patches.Rectangle((0.004,0.015), 0.012, 0.8, color='white', fill=True)
    #ax.add_patch(rect)
    #"""
    plt.text(0.00967, 10**-0.85, r'$n$ = {0:0.2f}'.format(m2_param[1]), fontsize=font_size, fontname="Times New Roman")
    #plt.text(0.00725, 10**-0.83, r'$R$ = %0.1f [min]' % ((1./(m2_param[2] / m2_param[0])**(-1./ m2_param[1]))/60.), fontsize=font_size, fontname="Times New Roman")
    #plt.text(0.007, 10**-1.09, r'$\beta$ = {0:0.3f}'.format(m2_param[4]), fontsize=25)
    plt.text(0.00975, 10**-1.1, r'$\beta$ = {0:0.1f} [min]'.format((1./np.exp(m2_param[4]))/60.), fontsize=font_size, fontname="Times New Roman")
    plt.text(0.01, 10**-1.35, r'$\delta$ = {0:0.3f}'.format(m2_param[5]), fontsize=font_size, fontname="Times New Roman")
    plt.text(0.01, 10**-1.6, r'$r$ = {0:0.3g}'.format(r_val[0]), fontsize=font_size, fontname="Times New Roman")
    #legend = ax.legend(loc='upper right', prop={'size':30}, labelspacing=0.35)
    legend = ax.legend(loc='lower left', prop={'size':font_size-4}, labelspacing=0.35)
    for label in legend.get_lines():
        label.set_linewidth(3.0)  # the legend line width
    plt.hold(False)
    #"""
    
    """
    plt.subplot(222)
    plt.loglog(freqs, avg_array)
    plt.vlines(1./164.,10**-3,10**-1)
    #plt.vlines(1./650.,10**-3,10**-1)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')
    plt.title('c) FFT Spectrum')
    plt.hold(False)
    """
    
    #--- Contour plot wavelet power spectrum
    plt3 = plt.subplot(223)
    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
    CS = plt.contourf(time, period, np.log2(power), len(levels))  #*** or use 'contour'
    im = plt.contourf(CS, levels=np.log2(levels))
    plt.xlabel('Time (seconds)')
    plt.ylabel('Period (seconds)')
    plt.title('b) Wavelet Power Spectrum (in base 2 logarithm)')
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
    
    wv_min = global_ws[0]
    fft_min = s[len(s)-1]
    adjust = wv_min/fft_min
    
    #--- Plot global wavelet spectrum
    plt4 = plt.subplot(224)
    ax = plt.gca()
    plt.hold(True)
    plt.ylim((10**-4.5,10**-0.5))
    plt.xlim((10**-4.,10**-1.5))
    plt.loglog(freqs, avg_array, 'k', label='FFT', linewidth=1.5)
    plt.loglog(1./(period*60.), global_ws/adjust, 'r', label='Wavelet', linewidth=1.5)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')
    plt.title('d) Global Wavelet Spectrum')
    plt.vlines(np.exp(m2_param[4]),10**-4.5,10**-0.5, 'b', linestyles='dashed', label='%0.2f minutes' % peak_fft, linewidth=1.5)
    #plt.vlines(1./650.,10**-1,10**2)
    #plt.xlim([0, 1.25 * np.max(global_ws)])
    # format y-scale
    #plt4.set_yscale('log', basey=2, subsy=None)
    #plt.ylim([np.min(period), np.max(period)])
    #ax = plt.gca().yaxis
    #ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #plt4.ticklabel_format(axis='y', style='plain')
    #plt4.invert_yaxis()
    legend = ax.legend(loc='lower left', prop={'size':font_size-4}, labelspacing=0.35)
    for label in legend.get_lines():
        label.set_linewidth(3.0)  # the legend line width
    plt.hold(False)
    plt.tight_layout()
    
    plt.show()
    plt.savefig('C:/Users/Brendan/Desktop/wavelet/wavelet_%ix%iy.jpeg' % (x,y),bbox_inches='tight')
    plt.close()