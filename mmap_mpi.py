# -*- coding: utf-8 -*-
"""
Created on Fri Feb 03 11:53:31 2017

@author: Brendan
"""

"""
### create memory-mapped array based on current array
"""

import numpy as np



# load original array  (I wrote down the first value in this array, to compare against later.)
original = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/spectra_20130530_1600_2300_2600i_2200_3000j_data_rebin4.npy')

# create memory-mapped array with similar datatype and shape to original array
mmap_arr = np.memmap('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/spectra_20130530_1600_2300_2600i_2200_3000j_data_rebin4_mmap.npy', dtype='float64', mode='w+', shape=(original.shape[0],original.shape[1],original.shape[2]))

# write data to memory-mapped array
mmap_arr[:] = original[:]

# flush memory changes to disk, then remove memory-mapped object
del mmap_arr



"""
### test that data was stored correctly 
"""

import numpy as np

# load memory-mapped array (have to specify shape - I guess record it somewhere, so don't have to load original to get its shape)
mmap_arr = np.memmap('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/spectra_20130530_1600_2300_2600i_2200_3000j_data_rebin4_mmap.npy', dtype='float64', mode='r', shape=(72,156,299))

# this value should match the first value in the original array
print mmap_arr[0][0][0] 



"""
### the codes I changed in the spec_fit script
"""

"""
### I changed these lines -- 
"""
cube = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_1600_2300_2600i_2200_3000j_data_rebin4_mmap.npy')  # works, (all load full array)


start = timer()

comm = MPI.COMM_WORLD 	# set up comms
rank = comm.Get_rank()	# Each processor gets its own "rank"


# Rank0 is the first processor. Use that to do the main chores
# possibly put cube load in rank == 0?, so only loaded once?
if rank == 0:
  size = MPI.COMM_WORLD.Get_size()		# How many processors do we have?
  #cube = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/spectra_20141025_304_-400_400i_-400_400j.npy')
  #cube = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130530/20130530_193_2300_2600i_2200_3000j_data_rebin1.npy') # works (1 array load)
  chunks = np.array_split(cube, size)		# Split the data based on no. of processors
else:
  chunks = None  # Prepare a variable on the other nodes

# Distribute the data to all nodes, defining the root node (rank=0) as distributor
subcube = comm.scatter(chunks, root=0)
ss = np.shape(subcube)  # Validation	
print "Processor", rank, "received an array with dimensions", ss  # Validation
print "Height = %i, Width = %i, Total pixels = %i" % (subcube.shape[0], subcube.shape[1], subcube.shape[0]*subcube.shape[1])



"""
### to these --
"""
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
	
start = timer()

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

# load memory-mapped array as read-only
cube = np.memmap('C:/Users/Brendan/Desktop/20130530_1600_2300_2600i_2200_3000j_data_rebin4_mmap.npy', dtype='float64', mode='r', shape=(72,156,299))

chunks = np.array_split(cube, size)  # Split the data based on no. of processors

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]

# verify each processor received subcube with correct dimensions
ss = np.shape(subcube)  # Validation	
print "Processor", rank, "received an array with dimensions", ss  # Validation
print "Height = %i, Width = %i, Total pixels = %i" % (subcube.shape[0], subcube.shape[1], subcube.shape[0]*subcube.shape[1])
