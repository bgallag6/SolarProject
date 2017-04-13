# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 19:56:25 2017

@author: Brendan
"""

from sunpy.net import vso
import astropy.units as u
import numpy as np
from mpi4py import MPI
import glob

def get_data_fill(arr_need, wavelength, directory):

    print "Please wait while request is being processed."

    client=vso.VSOClient()  # establish connection to VSO database
        
    # loop through the array of needed files, requesting them one at a time         
    for i in range(0,len(arr_need)):
        #qr=client.query(vso.attrs.Time(arr_need[i],arr_need[i]), vso.attrs.Instrument('aia'), vso.attrs.Wavelength(wavelength * u.AA, wavelength * u.AA))
        qr=client.query(vso.attrs.Time(arr_need[i],arr_need[i]), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
        #print qr
        #res=client.get(qr, path='%s/FITS/%s/%i/{file}.fits' % (directory, date, wavelength)).wait()
        res=client.get(qr, path='%s/FITS/%s/%i/{file}' % (directory, date, wavelength)).wait()  # use this from now on
        #print res
   
    flist = glob.glob('%s/FITS/%s/%i/*.fits' % (directory, date, wavelength))
    l = len(flist)
    arr_have = []
    l_fname = len(flist[0])
    for i in range(0,l):
        x = flist[i]
        h = int(x[(l_fname-28):(l_fname-26)])
        m = int(x[(l_fname-25):(l_fname-23)])
        s = int(x[(l_fname-22):(l_fname-20)])        
        t = ('%-11s''%02d'':''%02d'':''%02d' % (Y1,h,m,s))
        arr_have.append(t)
    arr_rem = list(set(arr_need) - set(arr_have))
    
    return arr_rem
 
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

import sys

directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])

#wavelength=193
#directory = 'F:/Users/Brendan/Desktop/SolarProject/FITS/20130530/193'    

arr_rem = np.load('%s/FITS/%s/%i/arr_need.npy' % (directory, date, wavelength))

Y1 = '%s/%s/%s' % (date[0:4],date[4:6],date[6:8])

counter = 0

for k in range(10):
    if len(arr_rem) > 0:
        if counter < 10:
            chunks = np.array_split(arr_rem, size)  # Split the data based on no. of processors
            sub_arr = comm.scatter(chunks, root=0)
            arr_rem_sub = get_data_fill(sub_arr, wavelength, directory)  # Do something with the array
            arr_rem_gat = comm.gather(arr_rem_sub, root=0)  # Gather all the results
            
            if rank == 0:
              arr_rem = sum(arr_rem_gat, [])
              print "Just finished pass %i" % counter
              print "Still need %i more files" % len(arr_rem)
              counter += 1

    