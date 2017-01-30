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

from timeit import default_timer as timer

import numpy as np
import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from matplotlib.widgets import  RectangleSelector
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
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
  params = np.zeros((7, SPECTRA.shape[0], SPECTRA.shape[1]))
  # M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], (len(freqs)+1)/2))  # would save storage / memory space
  M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], SPECTRA.shape[2]))

  # Uncertainties = np.zeros((6, SPECTRA.shape[0], SPECTRA.shape[1]))  # not using right now
  
  
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
        
        # populate array holding model fits
        M2_fit[l][m] = m2_fit
			
  return params, M2_fit
	

# load data
cube = np.load('C:/Users/Brendan/Desktop/SDO/20130530_1600_2300_2600i_2200_3000j_rebin4_spectra_mpi.npy')


start = timer()

comm = MPI.COMM_WORLD 	# set up comms
rank = comm.Get_rank()	# Each processor gets its own "rank"
#print "Hello World from process ", rank  # DEBUG/VALIDATE


# Rank0 is the first processor. Use that to do the main chores
if rank == 0:
  size = MPI.COMM_WORLD.Get_size()		# How many processors do we have?
  chunks = np.array_split(cube, size)		# Split the data based on no. of processors
else:
  chunks = None  # Prepare a variable on the other nodes

# Distribute the data to all nodes, defining the root node (rank=0) as distributor
subcube = comm.scatter(chunks, root=0)
ss = np.shape(subcube)  # Validation	
print "Processor", rank, "received an array with dimensions", ss  # Validation
print "Height = %i, Width = %i, Total pixels = %i" % (subcube.shape[0], subcube.shape[1], subcube.shape[0]*subcube.shape[1])
print "Estimated time remaining... "

params_T, M2_fit_T = spec_fit( subcube )		# Do something with the array
newData_p = comm.gather(params_T, root=0)	# Gather all the results
newData_m = comm.gather(M2_fit_T, root=0)	# Gather all the results

# Again, just have one node do the last bit
if rank == 0:
  #stack = np.vstack(newData) 	# stack the 2-d arrays together and we're done!
  stack_p = np.hstack(newData_p)
  #stack_m = np.vstack(newData_m)
  print stack_p.shape			# Verify we have a summed version of the input cube
  #print stack_m.shape			# Verify we have a summed version of the input cube
 
T_act = timer() - start
print "Program time = %i sec" % T_act     

np.save('C:/Users/Brendan/Desktop/SDO/20130530_1600_2300_2600i_2200_3000j_params_rebin4B', stack_p)
#np.save('C:/Users/Brendan/Desktop/SDO/M2_20130530_1600_2300_2600i_2200_3000j_data_rebin4_mpi_tool', stack_m)