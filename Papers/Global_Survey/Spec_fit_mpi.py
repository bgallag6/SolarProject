# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 12:15:32 2017

@author: Brendan
"""

"""
######################
# run with:
# $ mpiexec -n # python Spec_fit_mpi.py    (# = number of processors)
######################
"""

# maybe save results to text file - param bounds, region, fail-count... 

# 2/3:
# manually assign chunks to processors - to overcome 'Overflow error'

# 2/12
# changed ds from 0.1*s to new method

from timeit import default_timer as timer
import numpy as np
import scipy.signal
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
from scipy import fftpack
from astropy.convolution import convolve, Box1DKernel
from numpy.random import randn
from mpi4py import MPI
from scipy.stats.stats import pearsonr   

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)

# define combined-fitting function (Model M2)
def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)

def spec_fit( subcube ):
#def spec_fit( subcube, subcube_StdDev ):
    
  SPECTRA = subcube
  spectra_array = SPECTRA
  print SPECTRA.shape[0], SPECTRA.shape[1]
      
  # determine frequency values that FFT will evaluate
  num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used
  freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
  if wavelength == 1600 or wavelength == 1700:
      time_step = 24  # 24 second cadence for these wavelengths
      #time_step = 12  # Jacks dataset
  else:
      time_step = 12  # 12 second cadence for the others
      #time_step = 24  # for half-cadence test
  sample_freq = fftpack.fftfreq(freq_size, d=time_step)
  pidxs = np.where(sample_freq > 0)
  freqs = sample_freq[pidxs]


  # initialize arrays to hold parameter values
  params = np.zeros((9, SPECTRA.shape[0], SPECTRA.shape[1]))

  # Uncertainties = np.zeros((6, SPECTRA.shape[0], SPECTRA.shape[1]))  # not using right now
  
  start = timer()
  T1 = 0

    
  for l in range(SPECTRA.shape[0]):
  #for l in range(0,15):
    
    for m in range(SPECTRA.shape[1]):
    #for m in range(0,20):
        
                                        
        f = freqs  # frequencies
        s = spectra_array[l][m]  # fourier power
        
        
        # assign equal weights to all parts of the curve
        df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
        df2 = np.zeros_like(f)
        df2[0:len(df)] = df
        df2[len(df2)-1] = df2[len(df2)-2]
        ds = df2
        #ds = subcube_StdDev[l][m]
        
                                               
        ### fit data to models using SciPy's Levenberg-Marquart method
        
        try:
            # initial guesses for fitting parameters
            M1_low = [-0.002, 0.3, -0.01]
            M1_high = [0.002, 6., 0.01]
            nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
           
        
        except RuntimeError:
            #print("Error M1 - curve_fit failed - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        except ValueError:
            #print("Error M1 - inf/NaN - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass

      
        A, n, C = nlfit_l  # unpack fitting parameters

        ## fit data to combined power law plus gaussian component model
        #"""        
        try:                                 
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
            
            # change method to 'dogbox' and increase max number of function evaluations to 3000
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A,n,C,0.1,-5.55,0.425], bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
        
        except RuntimeError:
            #print("Error M2 - curve_fit failed - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        except ValueError:
            #print("Error M2 - inf/NaN - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        #"""
        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        #dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
        
        try:
            #if wavelength == 1600 or wavelength == 1700:
            #    M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]  # try constraining further, now that are specifying initial guess
            #    M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]
            #else:
            #    M2_low = [0., 0.3, 0., 0., -6.5, 0.05]  # try constraining further, now that are specifying initial guess
            #    M2_high = [0.002, 6., 0.01, 0.1, -4.6, 0.8]
            
            nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
            #nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000) # replaced #'s with arrays
       
        except RuntimeError:
            #print("Error M2 - curve_fit failed - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        except ValueError:
            #print("Error M2 - inf/NaN - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
        #dA22, dn22, dC22, dP22, dfp22, dfw22 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
        #m2_param = A22, n22, C22, P22, fp22, fw22  # could have used this for params array : = params[0:6,l-1,m-1]
        #uncertainties = dA22, dn22, dC22, dP22, dfp22, dfw22  # do we want to keep a global array of uncertainties?
                       
        # create model functions from fitted parameters
        m1_fit = PowerLaw(f, A, n, C)        
        #m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
        m2_fit2 = GaussPowerBase(f, A22,n22,C22,P22,fp22,fw22)      

        residsM1 = (s - m1_fit)
        chisqrM1 =  ((residsM1/ds)**2).sum()
        redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)  
                       
        #residsM2 = (s - m2_fit)
        #chisqrM2 = ((residsM2/ds)**2).sum()
        #redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
        
        residsM22 = (s - m2_fit2)
        chisqrM22 = ((residsM22/ds)**2).sum()
        redchisqrM22 = ((residsM22/ds)**2).sum()/float(f.size-6) 
        
        
        #f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))
        f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))
        
        #amp_scale = PowerLaw(np.exp(fp2), A2, n2, C2)  # to extract the gaussian-amplitude scaling factor
        amp_scale2 = PowerLaw(np.exp(fp22), A22, n22, C22)  # to extract the gaussian-amplitude scaling factor
        
        r_temp = pearsonr(m2_fit2, s)  # calculate r-value correlation coefficient
        r = r_temp[0]
        
        # populate array with parameters
        params[0][l][m] = A22
        params[1][l][m] = n22
        params[2][l][m] = C22
        params[3][l][m] = P22
        params[4][l][m] = fp22
        params[5][l][m] = fw22
        params[6][l][m] = f_test2
        params[7][l][m] = P22 / amp_scale2
        params[8][l][m] = r
        
        
    # estimate time remaining and print to screen  (looks to be much better - not sure why had above?)
    T = timer()
    T2 = T - T1
    if l == 0:
        T_init = T - start
        T_est = T_init*(SPECTRA.shape[0])  
        T_min, T_sec = divmod(T_est, 60)
        T_hr, T_min = divmod(T_min, 60)
        #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (l, SPECTRA.shape[0], T_est)
        print "Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (l, SPECTRA.shape[0], T_hr, T_min, T_sec)
    else:
        T_est2 = T2*((SPECTRA.shape[0])-l)
        T_min2, T_sec2 = divmod(T_est2, 60)
        T_hr2, T_min2 = divmod(T_min2, 60)
        #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (l, SPECTRA.shape[0], T_est2)
        print "Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (l, SPECTRA.shape[0], T_hr2, T_min2, T_sec2)
    T1 = T

  # print estimated and total program time to screen        
  print "Beginning Estimated time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec)
  T_act = timer() - start
  T_min3, T_sec3 = divmod(T_act, 60)
  T_hr3, T_min3 = divmod(T_min3, 60)
  print "Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3) 
			
  return params
	

comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
	
start = timer()

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

#directory = 'F:/Users/Brendan/Desktop/SolarProject'
#date = '20120923'
#wavelength = 94

import sys

# set variables from command line
directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])

# load memory-mapped array as read-only
cube_shape = np.load('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength))
cube = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))
#cube_StdDev = np.memmap('%s/DATA/Temp/%s/%i/uncertainties_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))
#cube = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20120923/171/20120923_171_-100_100i_-528_-132j_spectra.npy')

chunks = np.array_split(cube, size)  # Split the data based on no. of processors
#chunks_StdDev = np.array_split(cube_StdDev, size)  # Split the data based on no. of processors

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]
        #subcube_StdDev = chunks_StdDev[i]

# verify each processor received subcube with correct dimensions
ss = np.shape(subcube)  # Validation	
print "Processor", rank, "received an array with dimensions", ss  # Validation
print "Height = %i, Width = %i, Total pixels = %i" % (subcube.shape[0], subcube.shape[1], subcube.shape[0]*subcube.shape[1])

params_T = spec_fit( subcube )  # Do something with the array
#params_T = spec_fit( subcube, subcube_StdDev )  # Do something with the array
newData_p = comm.gather(params_T, root=0)  # Gather all the results

# Have one node stack the results
if rank == 0:
  stack_p = np.hstack(newData_p)
  print stack_p.shape  # Verify we have a summed version of the input cube
 

  T_final = timer() - start
  T_min_final, T_sec_final = divmod(T_final, 60)
  T_hr_final, T_min_final = divmod(T_min_final, 60)
  print "Total program time = %i:%.2i:%.2i" % (T_hr_final, T_min_final, T_sec_final)   
  print "Just finished region: %s %iA" % (date, wavelength)

  np.save('%s/DATA/Output/%s/%i/param' % (directory, date, wavelength), stack_p)