import numpy as np

# this function receives a 3D cube and sums over the z axis, returning the [x,y] sum array
def do_something( some_input ):
	xsize = np.shape( some_input )[0]
	ysize = np.shape( some_input )[1]
	sum_array = np.empty([xsize,ysize])
	#print sum_array.shape
	for i in range(0,xsize-1):
		for j in range(0,ysize-1):
			sum_array[i,j] = np.sum( some_input[i,j,:] )
			
	return sum_array
	
	
import numpy as np
from mpi4py import MPI

# generate some data
cube = np.load('spectra_20130815_193_1000_1600i_1950_2950j_rebin2.npy')

comm = MPI.COMM_WORLD 	# set up comms
rank = comm.Get_rank()	# Each processor gets its own "rank"
#print "Hello World from process ", rank		# DEBUG/VALIDATE


# Rank0 is the first processor. Use that to do the main chores
if rank == 0:
	size = MPI.COMM_WORLD.Get_size()		# How many processors do we have?
	chunks = np.array_split(cube, size)		# Split the data based on no. of processors
else:
	chunks = None			# Prepare a variable on the other nodes

# Distribute the data to all nodes, defining the root node (rank=0) as distributor
subcube = comm.scatter(chunks, root=0)
ss = np.shape(subcube)												# Validation	
print "Processor", rank, "received an array with dimensions", ss		# Validation

just_a_test = do_something( subcube )		# Do something with the array
newData = comm.gather(just_a_test,root=0)	# Gather all the results

# Again, just have one node do the last bit
if rank == 0:
	stack = np.vstack(newData) 	# stack the 2-d arrays together and we're done!
	print stack.shape			# Verify we have a summed version of the input cube
