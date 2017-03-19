# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 21:53:26 2017

@author: Brendan
"""

"""
############################
############################
# download data 
############################
############################
"""

# put in error handling if start date is after end date? also if data is entered incorrectly - not yyyy/mm/dd format

# generalize for 1600 - string bounds wont work - add 1

# maybe make first past be much quicker - like 10 min at a time

# for 20130815 193 - files are just ".fits" not ".fits.fits" - screws up call obviously

# 3/2: added extracting dates from query - removed manually creating the arr_all

from sunpy.net import vso
import astropy.units as u
import glob
import numpy as np

def get_data(wavelength, time_begin, time_end, path_name):
    """
    Downloads .fits image files from database. 
       
    wavelength : 
        Wavelength to be downloaded. (Integer)
       
    time_begin : 
        Beginning of time range in YYYY/MM/DD HH:MM:SS format. (String)
     
    time_end : 
        Ending of time range in YYYY/MM/DD HH:MM:SS format. (String)
                   
    path_name : 
        The directory to which the files should be downloaded. (String)
      
    Example:
    ::
        ss.get_data(wavelength=1600, time_begin='2016/09/23 00:00:00', time_end='2016/09/23 00:05:00',
           path_name='C:/Users/Brendan/Desktop/SDO_test')
    """
    
    print ""
    print "Please wait while request is being processed."

    client=vso.VSOClient()  # establish connection to VSO database
    
    # extract yyyy/mm/dd, hour, minute, and second values from start of time-range
    Y1 = str(time_begin[0:11])   
    H1 = int(time_begin[11:13])
    M1 = int(time_begin[14:16])
    S1 = int(time_begin[17:19])
    
    # extract yyyy/mm/dd, hour, minute, and second values from end of time-range
    Y2 = str(time_end[0:11])   
    H2 = int(time_end[11:13])
    M2 = int(time_end[14:16])
    S2 = int(time_end[17:19])
    
    # create time strings to pass to initial query request
    T1 = ('%s''%02d'':''%02d'':''%02d' % (Y1,H1,M1,S1))
    T2 = ('%s''%02d'':''%02d'':''%02d' % (Y2,H2,M2,S2))
    
    # query request to determine total number of files in time-range
    qr=client.query(vso.attrs.Time(T1,T2), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
    
    arr_all = [] 

    # extract all image times from query request to check against
    for i in range(len(qr)):
        sec = qr[i].time[1][12:14]
        minute = qr[i].time[1][10:12]
        hour = qr[i].time[1][8:10]
        day = qr[i].time[1][6:8]
        month = qr[i].time[1][4:6]
        year = qr[i].time[1][0:4]
        time = '%s/%s/%s' ' %s:%s:%s' % (year, month, day, hour, minute, sec)
        time = time.encode('utf8')
        arr_all = np.append(arr_all, time)
    
    num_files = len(qr)  # total number of files in time-range
    duration = (H2-H1)*3600 + (M2-M1)*60  # total length of time-range
    spacing = duration / num_files  # determine average image cadence
    if spacing > 10 and spacing < 14:
        cadence = 12
    elif spacing > 22 and spacing < 26:
        cadence = 24
    
    print ""    
    print "This request will download the dataset of wavelength %d" % wavelength
    print "in the timerange from %s to %s." % (T1,T2)
    print "Image files will be saved in the directory: %s" % path_name
    print ""
    print "This will download %d files at %d-second cadence." % (num_files,cadence)
    print ""
    raw_input("Press [Enter] to continue...")
    
    # set starting values of iteration to specified start time
    s1 = S1
    s2 = s1
    m1 = M1
    m2 = m1 + 1
    #m2 = m1 + mins  # if want to try for faster 
    h1 = H1
    h2 = h1
    
    # loop through entire time range, 1 minute at a time (so as not to overload the server with requests)
    for i in range((duration/60)-1):
    #for i in range(0,(duration/(60*mins))-1):  # try for quicker download using custom minutes-per-query
        m1 += 1
        m2 += 1
        #m1 += mins
        #m2 += mins
        if m1 >= 60:
            m1 -= 60
            h1 += 1
        if m2 >= 60:
            m2 -= 60
            h2 += 1
        t1 = ('%s''%02d'':''%02d'':''%02d' % (Y1,h1,m1,s1))
        print t1
        t2 = ('%s''%02d'':''%02d'':''%02d' % (Y1,h2,m2,s2))
        print t2
        
        qr=client.query(vso.attrs.Time(t1,t2), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
        #print qr
        res=client.get(qr, path='%s/{file}.fits' % path_name).wait()  # leave the "{file}.fits" part alone  <-- ? 
    
    
    ## after initial passthrough, determine files downloaded and those that still need to be    
        
    # get list of files downloaded    
    flist = glob.glob('%s/*.fits' % path_name)
    
    l = len(flist)
     
    l_fname = len(flist[0])
    
    arr_have = []
    
    # create searchable array of images that have already been downloaded
    for i in range(l):
        x = flist[i]
        h = int(x[(l_fname-33):(l_fname-31)])
        m = int(x[(l_fname-30):(l_fname-28)])
        s = int(x[(l_fname-27):(l_fname-25)])
        t = ('%s''%02d'':''%02d'':''%02d' % (Y1,h,m,s))
        arr_have.append(t)
    #print arr_have
    
    # compare array_all to array_have to determine array_need
    arr_need = []
    for i in range(num_files):  # changed from range(0,l)
        z = arr_all[i] in arr_have    
        if z == False:
            arr_need.append(arr_all[i])
    #print arr_need
    
    print ""
    print "After the initial pass, still need %d files." % len(arr_need)
    
    # loop through the array of needed files, requesting them one at a time        
    for i in range(len(arr_need)):
        qr=client.query(vso.attrs.Time(arr_need[i],arr_need[i]), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
        res=client.get(qr, path='%s/{file}.fits' % path_name).wait()  # .wait() -- wait until file is downloaded before issuing another request
  
      

"""
############################
############################
# re-download data 
############################
############################
"""        
        
### maybe put as an if statement in the original - or specify argument in original call
# like (first / new   vs re / fill)
        
# prof weigel mentioned that could instead use push/pop to move first element of array to bottom
# so that program doesn't get stuck on one file, also that it would be better to use the ISO module for date/time

        
from sunpy.net import vso
import astropy.units as u
import glob
import numpy as np

def get_data_fill(wavelength, time_begin, time_end, path_name):
    """
    Downloads .fits image files from database. 
    
    wavelength : 
        Wavelength to be downloaded. (Integer)
       
    time_begin : 
        Beginning of time range in YYYY/MM/DD HH:MM:SS format. (String)
     
    time_end : 
        Ending of time range in YYYY/MM/DD HH:MM:SS format. (String)
                   
    path_name : 
        The directory to which the files should be downloaded. (String)
      
    Example:
    ::
        ss.get_data_fill(wavelength=171, time_begin='2013/08/15 00:00:00',
           time_end='2013/08/15 12:00:00', path_name='F:/SDO/data/20130815/171')
    """
    
    print ""
    print "Please wait while request is being processed."

    client=vso.VSOClient()  # establish connection to VSO database
    
    # extract yyyy/mm/dd, hour, minute, and second values from start of time-range
    Y1 = str(time_begin[0:11])   
    H1 = int(time_begin[11:13])
    M1 = int(time_begin[14:16])
    S1 = int(time_begin[17:19])
    
    # extract yyyy/mm/dd, hour, minute, and second values from end of time-range
    Y2 = str(time_end[0:11])   
    H2 = int(time_end[11:13])
    M2 = int(time_end[14:16])
    S2 = int(time_end[17:19])
    
    # create time strings to pass to initial query request
    T1 = ('%s''%02d'':''%02d'':''%02d' % (Y1,H1,M1,S1))
    T2 = ('%s''%02d'':''%02d'':''%02d' % (Y2,H2,M2,S2))
    
    # query request to determine total number of files in time-range
    qr=client.query(vso.attrs.Time(T1,T2), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
    
    arr_all = [] 

    for i in range(len(qr)):
        sec = qr[i].time[1][12:14]
        minute = qr[i].time[1][10:12]
        hour = qr[i].time[1][8:10]
        day = qr[i].time[1][6:8]
        month = qr[i].time[1][4:6]
        year = qr[i].time[1][0:4]
        time = '%s/%s/%s' ' %s:%s:%s' % (year, month, day, hour, minute, sec)
        time = time.encode('utf8')
        arr_all = np.append(arr_all, time)
    
    num_files = len(qr)
    print num_files
    #cadence = cadence  # set cadence to specified value *took out because extracting times from query
        
    flist = glob.glob('%s/*.fits' % path_name)
    
    l = len(flist)
     
    l_fname = len(flist[0])
   
    # find first file after 00:00:00 - set time base to that
    
    # loop through flist couple times until get all.
    
    
    # create searchable array of images that have already been downloaded
        
    adj = 5  # adjust 5 characters for 20130815 193 dataset (doesn't have extra '.fits')
    #adj = 0  # for all other datasets
    
    arr_have = []
    for i in range(0,l):
        x = flist[i]
        h = int(x[(l_fname-33+adj):(l_fname-31+adj)])
        m = int(x[(l_fname-30+adj):(l_fname-28+adj)])
        s = int(x[(l_fname-27+adj):(l_fname-25+adj)])        
        t = ('%s''%02d'':''%02d'':''%02d' % (Y1,h,m,s))
        arr_have.append(t)
    #print arr_have
    print len(arr_have)    
    
    # compare array_all to array_have to determine array_need
    arr_need = []   
    for i in range(0,num_files):  # changed from range(0,l)
        z = arr_all[i] in arr_have    
        if z == False:
            arr_need.append(arr_all[i])
    print arr_need
    print len(arr_need)    
    
    print ""
    print "After the initial pass, still need %d files." % len(arr_need)
    
    # loop through the array of needed files, requesting them one at a time         
    for i in range(0,len(arr_need)):
        qr=client.query(vso.attrs.Time(arr_need[i],arr_need[i]), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
        print qr
        #res=client.get(qr, path='%s/{file}.fits' % path_name).wait()
        res=client.get(qr, path='%s/{file}' % path_name).wait()  # only for 193  (maybe for all going forward?)
        #print res
        


"""
############################
############################
# Curve Fitting
############################
############################
"""    

# should revise the variable naming ('gp' --> 'M2', all the '2' variables?)

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import h5py
from scipy import fftpack  # doesnt work in module when called here???
import matplotlib.pylab as plt
from astropy.convolution import convolve, Box1DKernel
from matplotlib import cm
from numpy.random import randn
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.ticker import LogFormatterMathtext
from timeit import default_timer as timer
from scipy.stats import f


def spec_fit(spectra_array):
    """
    Calculates pixel-box averages and fits spectra.
    
    spectra_array : 
        3D array of FFT-segment and 3x3-pixel-box averaged region  (array)
        
    Returns : 
        Parameter and M2-fit arrays
      
    Example:
    ::
        ss.spec_fit(spectra_array = SPECTRA)
    """

    from scipy import fftpack    
    
    # define Power-Law-fitting function (Model M1)
    def PowerLaw(f, A, n, C):
        return A*f**-n + C
        
    # define Gaussian-fitting function
    def Gauss(f, P, fp, fw):
        return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)
    
    # define combined-fitting function (Model M2)
    def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
        return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)
        
    
    ## load in array of segment-averaged pixel FFTs
    SPECTRA = spectra_array
    #spectra_array = spectra_array.astype(np.float32)  # possibly use this?
    
    print "The region size is %ii x %ij" % (SPECTRA.shape[0], SPECTRA.shape[1])
    print "%i frequencies were evaluated in the FFT" % SPECTRA.shape[2] 
    
    num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used
        
    # determine frequency values that FFT will evaluate
    freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
    time_step = 12  # add as argument, or leave in as constant?
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    #freqs = freqs.astype(np.float32)  # possibly use this?
    
    
    # initialize arrays to hold parameter values, also each pixel's combined model fit - for tool
    #diffM1M2 = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1]))  # dont really use - get rid of?
    params = np.zeros((8, SPECTRA.shape[0], SPECTRA.shape[1]))
    #M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], (len(freqs)+1)/2))  # would save storage / memory space
    M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], SPECTRA.shape[2]))
    
    #Uncertainties = np.zeros((6, SPECTRA.shape[0], SPECTRA.shape[1]))
    
    start = timer()
    T1 = 0
    
    
    ### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
    ## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)
    
    for l in range(0,SPECTRA.shape[0]):
    #for l in range(145,146):
        
        for m in range(0,SPECTRA.shape[1]):
        #for m in range(185,195):
            
                                            
            f = freqs  # frequencies
            s = spectra_array[l][m]  # fourier power
            
            #ds = (1./f**2.2)/1000
            #ds = s*0.1  # set the error / variance estimate to a constant percentage of the spectra power-values
            
            # assign equal weights to all parts of the curve
            df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
            df2 = np.zeros_like(f)
            df2[0:len(df)] = df
            df2[len(df2)-1] = df2[len(df2)-2]
            ds = df2
            
            # create points to fit model with final parameters 
            #f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],(len(freqs)+1)/2)  # would save storage / memory space?
            #f_fit = freqs       
            
                                                   
            ### fit data to models using SciPy's Levenberg-Marquart method
            
            try:
                # initial guesses for fitting parameters
                M1_low = [-0.002, 0.3, -0.01]
                M1_high = [0.002, 6., 0.01]
                nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
               
            
            except RuntimeError:
                print("Error M1 - curve_fit failed - %i, %i" % (l,m))
            
            except ValueError:
                print("Error M1 - inf/NaN - %i, %i" % (l,m))
    
          
            A, n, C = nlfit_l  # unpack fitting parameters
            
            # unpack uncertainties in fitting parameters from diagonal of covariance matrix
            dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
            
            
            ## possibly use previous pixel's parameters as initial guesses for current pixel (issues creating wierd banding in images)
                            
            """
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            #M2_low = [-0.1, 0.1, -0.1, 0.00001, -6.5, 0.05]  # test on 193 - coronal hole
            M2_high = [0.002, 4., 0.01, 0.2, -4.6, 0.8]
            if m > 0:
                P0 = [params[0][l][m-1], params[1][l][m-1], params[2][l][m-1], params[3][l][m-1], params[4][l][m-1], params[5][l][m-1]]
                #P0 = [0.0, 1.02, 0.001, 0.001, -4.68, 0.79]
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
                    
            elif m == 0 and l > 0:
                P0 = [params[0][l-1][m], params[1][l-1][m], params[2][l-1][m], params[3][l-1][m], params[4][l-1][m], params[5][l-1][m]]
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
                    
            else:        
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0=[nlfit_l[0], nlfit_l[1], nlfit_l[2], nlfit_g[0], nlfit_g[1], nlfit_g[2]], sigma=ds)
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
            """     
            
            ## fit data to combined power law plus gaussian component model
            #"""        
            try:                                 
                M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
                M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
                #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
                
                # change method to 'dogbox' and increase max number of function evaluations to 3000
                nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
                
            except RuntimeError:
                print("Error M2 - curve_fit failed - %i, %i" % (l,m))
            
            except ValueError:
                print("Error M2 - inf/NaN - %i, %i" % (l,m))
            #"""
            
            A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
            
            # unpack uncertainties in fitting parameters from diagonal of covariance matrix
            dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
            
            m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]
            uncertainties = dA2, dn2, dC2, dP2, dfp2, dfw2  # do we want to keep a global array of uncertainties?  
            
            
            #uncertainties_arr = [dA2, dn2, dC2, dP2, dfp2, dfw2]  # not sure if want to keep these
            #Uncertainties[:, l, m] = uncertainties_arr
            
            
            # create model functions from fitted parameters
            #m1_fit = PowerLaw(f_fit, A, n, C)
            m1_fit = PowerLaw(f, A, n, C)
            amp_scale = PowerLaw(np.exp(fp2), A, n, C)  # to extract the gaussian-amplitude scaling factor
            #m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
            m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
            #s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
            #m2P_fit = PowerLaw(f_fit, A2, n2, C2)  # only need if plotting
            #m2G_fit = Gauss(f_fit, P2, fp2, fw2)  # only need if plotting
            
            #diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
            #diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
                                   
            
            #residsM2 = (s - s_fit_gp_full)
            residsM2 = (s - m2_fit)
            chisqrM2 = ((residsM2/ds)**2).sum()
            redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
            
            residsM1 = (s - m1_fit)
            chisqrM1 =  ((residsM1/ds)**2).sum()
            redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       
            
            f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))
            
            # populate array with parameters
            params[0][l][m] = A2
            params[1][l][m] = n2
            params[2][l][m] = C2
            params[3][l][m] = P2
            params[4][l][m] = fp2
            params[5][l][m] = fw2
            #params[6][l][m] = redchisqrM2
            params[6][l][m] = f_test
            params[7][l][m] = P2 / amp_scale
            
            
            # populate array holding model fits
            M2_fit[l][m] = m2_fit
            
            
            
            # Plot models + display combined-model parameters + uncertainties
            """
            fig = plt.figure(figsize=(20,15))
            plt.title('SDO AIA 304.0 Angstrom 20120923 - 3 Segments, 3x3 Pixel-Box Averaging 598 interp', y = 1.01, fontsize=25)
            plt.ylim((10**-8,10**1))
            plt.xlim((10**-5,10**-1))
            plt.loglog(f,s,'k')
            plt.loglog(f_fit, m1_fit, label='Power Law - M1')
            plt.loglog(f_fit, m2P_fit, 'g', label='Power Law - M2')
            plt.loglog(f_fit, m2G_fit, 'g--', label='Gaussian - M2')
            plt.loglog(f_fit, m2_fit, 'r', label='Combined - M2')
            plt.xlabel('Frequency (Hz)', fontsize=20, labelpad=10)
            plt.ylabel('Power', fontsize=20, labelpad=10)
            plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
            plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
            plt.vlines((1.0/24.),10**-8,10**1, linestyles='solid', label='24 seconds')
            plt.text(0.01, 10**0., 'A = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[0], uncertainties[0]), fontsize=15)
            plt.text(0.01, 10**-0.3, 'n = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[1], uncertainties[1]), fontsize=15)
            plt.text(0.01, 10**-0.6, 'C = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[2], uncertainties[2]), fontsize=15)
            plt.text(0.01, 10**-0.9, 'P = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[3], uncertainties[3]), fontsize=15)
            plt.text(0.01, 10**-1.2, 'fp = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[4], uncertainties[4]), fontsize=15)
            plt.text(0.01, 10**-1.5, 'fw = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[5], uncertainties[5]), fontsize=15)
            plt.text(0.01, 10**-1.8, 'f_test = {0:0.3f}'.format(f_test), fontsize=15)
            plt.legend(loc='upper left', prop={'size':15})
            plt.show()
            #plt.savefig('C:/Users/Brendan/Desktop/PHYS 326/dogbox_test1/20130530_193A_3x3_6seg_%ii_%ij.jpeg' % (l,m))
            #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
            #plt.close()
            """
        
        # estimate time remaining and print to screen  (looks to be much better - not sure why had above?)
        T = timer()
        T2 = T - T1
        if l == 0:
            T_init = T - start
            T_est = T_init*(SPECTRA.shape[0])  
            print "Currently on row %i of %i, estimated time remaining: %i seconds" % (l, SPECTRA.shape[0], T_est)
        else:
            T_est2 = T2*((SPECTRA.shape[0])-l)
            print "Currently on row %i of %i, estimated time remaining: %i seconds" % (l, SPECTRA.shape[0], T_est2)
        T1 = T
    
    # print estimated and total program time to screen        
    print "Beginning Estimated time = %i sec" % T_est
    T_act = timer() - start
    print "Actual total time = %i sec" % T_act           
    
    return params, M2_fit
    
    


