# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 04:01:07 2017

@author: Brendan
"""

from sunpy.net import vso
import astropy.units as u
import numpy as np
from mpi4py import MPI

def get_data_fill(arr_need, wavelength, directory):

    print "Please wait while request is being processed."

    client=vso.VSOClient()  # establish connection to VSO database
        
    # loop through the array of needed files, requesting them one at a time         
    for i in range(0,len(arr_need)):
        qr=client.query(vso.attrs.Time(arr_need[i],arr_need[i]), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
        print qr
        #res=client.get(qr, path='%s/FITS/%s/%i/{file}.fits' % (directory, date, wavelength)).wait()
        res=client.get(qr, path='%s/FITS/%s/%i/{file}' % (directory, date, wavelength)).wait()  # use this from now on
        #print res
 
 
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

import sys

directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])

#wavelength=193
#directory = 'F:/Users/Brendan/Desktop/SolarProject/FITS/20130530/193'    

arr_need = np.load('%s/FITS/%s/%i/arr_need.npy' % (directory, date, wavelength))

chunks = np.array_split(arr_need, size)  # Split the data based on no. of processors

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        sub_arr = chunks[i]


get_data_fill(sub_arr, wavelength, directory)  # Do something with the array