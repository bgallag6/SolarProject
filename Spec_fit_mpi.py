# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 12:15:32 2017

@author: Brendan
"""

"""
######################
# run with:
# $ mpiexec -n 8 python Spec_fit_mpi.py    (8 = number of processors)
######################
"""

# dont absolutely need M2_fit
# can get rid of M2_gauss, M2_powerlaw, s_fit_gp

# param + m2_fit both generated perfectly

# 1/29:
# took out f_fit, replaced with f, since were the same

# maybe print param bounds?

# maybe save results to text file - param bounds, region, fail-count... 

# commenting-out the M2_fits.  for large regions they become way to big.  
# include back when need them?

# 2/3:
# manually assign chunks to processors - to overcome 'Overflow error'

from timeit import default_timer as timer

import numpy as np
import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
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
from mpi4py import MPI

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

# this function receives a 3D cube and sums over the z axis, returning the [x,y] sum array
def spec_fit( subcube ):
    
  SPECTRA = subcube
  spectra_array = SPECTRA
  print SPECTRA.shape[0], SPECTRA.shape[1]
  num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used
    
  # determine frequency values that FFT will evaluate
  freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
  time_step = 12  # add as argument, or leave in as constant?
  sample_freq = fftpack.fftfreq(freq_size, d=time_step)
  pidxs = np.where(sample_freq > 0)
  freqs = sample_freq[pidxs]
  #freqs = freqs.astype(np.float32)  # possibly use this?


  # initialize arrays to hold parameter values, also each pixel's combined model fit - for tool
  # diffM1M2 = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1]))  # not using right now
  params = np.zeros((8, SPECTRA.shape[0], SPECTRA.shape[1]))
  # M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], (len(freqs)+1)/2))  # would save storage / memory space
  #M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], SPECTRA.shape[2]))

  # Uncertainties = np.zeros((6, SPECTRA.shape[0], SPECTRA.shape[1]))  # not using right now
  
  start = timer()
  T1 = 0
  
  ### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
  ## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)
    
  for l in range(0,SPECTRA.shape[0]):
  #for l in range(0,15):
    
    for m in range(0,SPECTRA.shape[1]):
    #for m in range(0,20):
        
                                        
        f = freqs  # frequencies
        s = spectra_array[l][m]  # fourier power
        
        #ds = (1./f**2.2)/1000
        ds = s*0.1  # set the error / variance estimate to a constant percentage of the spectra power-values
        
        # create points to fit model with final parameters 
        #f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],(len(freqs)+1)/2)  # would save storage / memory space?
        #f_fit = freqs       
        
                                               
        ### fit data to models using SciPy's Levenberg-Marquart method
        
        try:
            # initial guesses for fitting parameters
            M1_low = [-0.002, 0.3, -0.01]
            M1_high = [0.002, 4., 0.01]
            nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
           
        
        except RuntimeError:
            #print("Error M1 - curve_fit failed - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass
        
        except ValueError:
            #print("Error M1 - inf/NaN - %i, %i" % (l,m))  # turn off because would print too many to terminal
            pass

      
        A, n, C = nlfit_l  # unpack fitting parameters
        
        # unpack uncertainties in fitting parameters from diagonal of covariance matrix
        dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]

        ## fit data to combined power law plus gaussian component model
        #"""        
        try:                                 
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            M2_high = [0.002, 4., 0.01, 0.2, -4.6, 0.8]
            #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
            
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
        dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
        
                   
        # create model functions from fitted parameters
        #m1_fit = PowerLaw(f_fit, A, n, C)
        m1_fit = PowerLaw(f, A, n, C)
        amp_scale = PowerLaw(np.exp(fp2), A, n, C)  # to extract the gaussian-amplitude scaling factor
        #m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
        m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
        #s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
        # m2P_fit = PowerLaw(f_fit, A2, n2, C2)  # only need if plotting
        # m2G_fit = Gauss(f_fit, P2, fp2, fw2)  # only need if plotting 
        
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
        #params[8][l][m] = 
        
        # populate array holding model fits
        #M2_fit[l][m] = m2_fit
        
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
			
  #return params, M2_fit
  return params
	

comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
	
start = timer()

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

# load memory-mapped array as read-only
cube = np.memmap('F:/Users/Brendan/Desktop/SolarProject/data/20130626/211/20130626_211_-500_500i_-500_500j_spectra_mmap.npy', dtype='float64', mode='r', shape=(1658,1481,299))

chunks = np.array_split(cube, size)  # Split the data based on no. of processors

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]

# verify each processor received subcube with correct dimensions
ss = np.shape(subcube)  # Validation	
print "Processor", rank, "received an array with dimensions", ss  # Validation
print "Height = %i, Width = %i, Total pixels = %i" % (subcube.shape[0], subcube.shape[1], subcube.shape[0]*subcube.shape[1])

#params_T, M2_fit_T = spec_fit( subcube )		# Do something with the array
params_T = spec_fit( subcube )		# Do something with the array
newData_p = comm.gather(params_T, root=0)	# Gather all the results
#newData_m = comm.gather(M2_fit_T, root=0)	# Gather all the results

# Again, just have one node do the last bit
if rank == 0:
  #stack = np.vstack(newData) 	# stack the 2-d arrays together and we're done!
  stack_p = np.hstack(newData_p)
  #stack_m = np.vstack(newData_m)
  print stack_p.shape			# Verify we have a summed version of the input cube
  #print stack_m.shape			# Verify we have a summed version of the input cube
 
T_act = timer() - start
print "Program time = %i sec" % T_act     

#np.save('/media/brendan/My Passport/Users/Brendan/Desktop/SolarProject/spectra_20130815_193_1000_1600i_1950_2950j_rebin2_params_mpi', stack_p)
np.save('F:/Users/Brendan/Desktop/SolarProject/data/20130626/211/20130626_211_-500_500i_-500_500j_param_mmap', stack_p)