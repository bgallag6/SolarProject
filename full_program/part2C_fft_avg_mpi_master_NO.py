# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 15:31:24 2018

@author: Brendan
"""

"""
######################
# run with:
# $ mpiexec -n # python part3B_spec_fit_mpi.py    (# = number of processors)
######################
"""


from timeit import default_timer as timer
import numpy as np
import scipy.signal
from scipy import fftpack
from mpi4py import MPI
import yaml

DIETAG = 9999 

class Work():
    def __init__(self, work_items):
        self.work_items = work_items[:] 
 
    def get_next_item(self):
        if len(self.work_items) == 0:
            return None
        return self.work_items.pop()                 

def fft_avg(subcube):
        
    pixmed=np.empty(subcube.shape[0])  # Initialize array to hold median pixel values
    spectra_seg = np.zeros((subcube.shape[1],len(freqs)))
        
    for col in range(subcube.shape[1]):       
        
            pixmed = subcube[:,col] / exposure  # extract timeseries + normalize by exposure time   
        
            v_interp = np.interp(t_interp,time,pixmed)  # interpolate pixel-intensity values onto specified time grid
            
            data = v_interp
            
            avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
    
            data = data[:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
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
            
            avg_array /= n_segments  # take the average of the segments
            
            avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
                       
            spectra_seg[col] = avg_array  # construct 3D array with averaged FFTs from each pixel
                
    return spectra_seg

 
def master():
    all_data = np.zeros((cube_shape[1],cube_shape[2],len(freqs)))
    cube = np.memmap('%s/DATA/Temp/%s/%i/derotated_mmap.npy' % (directory, date, wavelength), dtype='int16', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))
    workload = list(np.arange(cube_shape[1]-1,-1,-1))
    print(len(workload))
    
    size = MPI.COMM_WORLD.Get_size()
    current_work = Work(workload) 
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    
    start = timer()
    T1 = 0
    counter = 0
    counter0 = 0
    
    for i in range(1, size): 
        rnext = current_work.get_next_item() 
        if rnext == None: break
        print(cube[:,rnext].shape[0],cube[:,rnext].shape[1])
        comm.send(obj=cube[:,rnext], dest=i, tag=rnext)
 
    while 1:
        rnext = current_work.get_next_item()
        if rnext == None: break
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        tg = status.Get_tag()
        all_data[tg,:] = data

        comm.send(obj=cube[:,rnext], dest=status.Get_source(), tag=rnext)
        
        # estimate time remaining and print to screen
        if status.Get_source() == 1:
            T = timer()
            T2 = T - T1
            if counter0 == 0:
                T_init = T - start
                T_est = T_init*(cube_shape[1])/(size-1)  
                T_min, T_sec = divmod(T_est, 60)
                T_hr, T_min = divmod(T_min, 60)
                print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (counter+1, cube_shape[1], T_hr, T_min, T_sec), flush=True)
            else:
                T_est2 = T2*(cube_shape[1]-counter)/(size-1)
                T_min2, T_sec2 = divmod(T_est2, 60)
                T_hr2, T_min2 = divmod(T_min2, 60)
                print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (counter+1, cube_shape[1], T_hr2, T_min2, T_sec2), flush=True)
            T1 = T
            counter0 = 1
        counter +=1
 
    for i in range(1,size):
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        tg = status.Get_tag()
        all_data[tg,:] = data
    
    for i in range(1,size):
        comm.send(obj=None, dest=i, tag=DIETAG)
        
    # print estimated and total program time to screen        
    print("Beginning Estimated time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec), flush=True)
    T_act = timer() - start
    T_min3, T_sec3 = divmod(T_act, 60)
    T_hr3, T_min3 = divmod(T_min3, 60)
    print("Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3), flush=True) 
     
    return all_data
        
    
def slave():
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    while 1:
        data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        if status.Get_tag() == DIETAG: break
        p = fft_avg(data)
        comm.send(obj=p, dest=0, tag=status.Get_tag())
  

"""
# Setup MPI and load datacube, time, and exposure arrays
"""  

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)
rank = MPI.COMM_WORLD.Get_rank()  # Each processor gets its own "rank"

stream = open('specFit_config.yaml', 'r')
cfg = yaml.load(stream)

directory = cfg['temp_dir']
date = cfg['date']
wavelength = cfg['wavelength']
mmap_spectra = cfg['mmap_spectra']
save_temp = cfg['save_temp']
n_segments = cfg['num_segments']  # break data into # segments of equal length

cube_shape = np.load('%s/DATA/Temp/%s/%i/derotated_mmap_shape.npy' % (directory, date, wavelength))

time = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))
exposure = np.load('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength))

# determine frequency values that FFT will evaluate
if wavelength in [1600,1700]:
    time_step = 24  # add as argument in function call, or leave in as constant?
else:
    time_step = 12

t_interp = np.linspace(0, time[len(time)-1], (time[len(time)-1]//time_step)+1)  # interpolate onto default-cadence time-grid
    
#n_segments = num_seg
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) // n_segments

sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]


if rank == 0:
    all_dat = master()
    
    """
    ### 3x3 Averaging
    """
    temp = np.zeros((9,len(freqs)))  # maybe have 3x3 to be generalized   
    spectra_array = np.zeros((cube_shape.shape[1]-2, cube_shape.shape[2]-2, len(freqs)))
    spectra_StdDev = np.zeros((cube_shape.shape[1]-2, cube_shape.shape[2]-2, len(freqs)))
    
    ### calculate 3x3 pixel-box average, start at 1 and end 1 before to deal with edges
    
    for l in range(1,cube_shape.shape[1]-1):
    
        for m in range(1,cube_shape.shape[2]-1):
            
            temp[0] = all_dat[l-1][m-1]
            temp[1] = all_dat[l-1][m]
            temp[2] = all_dat[l-1][m+1]
            temp[3] = all_dat[l][m-1]
            temp[4] = all_dat[l][m]
            temp[5] = all_dat[l][m+1]
            temp[6] = all_dat[l+1][m-1]
            temp[7] = all_dat[l+1][m]
            temp[8] = all_dat[l+1][m+1]
    
            p_avg = np.average(temp, axis=0)
            
            spectra_array[l-1][m-1] = p_avg
            spectra_StdDev[l-1][m-1] = np.std(temp, axis=0)    
            
    print("Saving parameter file...", flush=True)
    
    if mmap_spectra == "y":
        orig_shape = np.array([spectra_array.shape[0], spectra_array.shape[1], spectra_array.shape[2]])
        
        # create memory-mapped array with similar datatype and shape to original array
        mmap_arr = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='%s' % spectra_array.dtype, mode='w+', shape=tuple(orig_shape))
        mmap_StdDev = np.memmap('%s/DATA/Temp/%s/%i/3x3_stddev_mmap.npy' % (directory, date, wavelength), dtype='%s' % spectra_array.dtype, mode='w+', shape=tuple(orig_shape))
        
        # write data to memory-mapped array
        mmap_arr[:] = spectra_array[:]
        mmap_StdDev[:] = spectra_StdDev[:]
        
        # save memory-mapped array dimensions to use when loading
        np.save('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength), orig_shape)
    
        # save original array if specified
        if save_temp == "y":
            np.save('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength), spectra_array)
            np.save('%s/DATA/Temp/%s/%i/3x3_stddev.npy' % (directory, date, wavelength), spectra_StdDev)
        
        # flush memory changes to disk, then remove memory-mapped object and original array
        del mmap_arr
        del mmap_StdDev
        del spectra_array
        del spectra_StdDev        
        
    else:
        np.save('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength), spectra_array)
        np.save('%s/DATA/Temp/%s/%i/3x3_stddev.npy' % (directory, date, wavelength), spectra_StdDev)
      
    print("Just finished region: %s %iA" % (date, wavelength), flush=True)
else:
    slave() 