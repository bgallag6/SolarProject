# -*- coding: utf-8 -*-
"""
Created on Mon Oct 02 19:04:02 2017

@author: Brendan
"""

import numpy as np
from mpi4py import MPI
import urllib.request
#import urllib  # Python 2
#import urllib2  # Python 2
import os

#"""
def get_data_fill(arr_need, arr_rename, directory):
    
    if rank == 0:
        print("Please wait while request is being processed.")    
        sub_arr_len = len(arr_need)
    
    counter = 0
        
    # loop through the array of needed files, requesting them one at a time         
    for i in range(len(arr_need)):
        
        if url_want.find("hmi.M") != -1:
            wavelength = 'magnetogram'

        elif url_want.find("hmi.I") != -1:
            wavelength = 'continuum'
        
        elif url_want.find("aia.l") != -1:
            wavelength_temp = arr_need[i][-20:-16]
            #wavelength = arr_need[i][-19:-16]
            if wavelength_temp[0] == '.':
                wavelength = wavelength_temp[1:4]
            else:
                wavelength = wavelength_temp
            #print wavelength
                
        wavelengths = ['magnetogram', 'continuum', '1700', '1600', '304', '171', '193', '211']
        
        if wavelength in wavelengths:
            
            fn=("%s/FITS/%s/%s/%s" % (directory, date, wavelength, arr_rename[i]))

            if not os.path.isfile(fn):
                urllib.request.urlretrieve("%s" % arr_need[i], fn) 
            
            counter += 1
        
        if rank == 0:
            print("Probably downloading file %i/%i" % (counter*size, sub_arr_len*size))
    return counter


              
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

import sys

#jsoc_url = sys.argv[1]
#directory = sys.argv[2]
#date = sys.argv[3]
directory = sys.argv[1]
date = sys.argv[2]
wavelength = sys.argv[3]

r_url = np.load('%s/FITS/%s/%s_request_url.npy' % (directory, date, wavelength))  # only use with JSOC_request_url.py - otherwise switch commented stuff back
arr_need = []
arr_rename = []

#page=urllib2.urlopen(jsoc_url)
#page=urllib2.urlopen(str(r_url))  # Python 2
#data=page.read().split("<td><A HREF=")  # Python 2
page=urllib.request.urlopen(str(r_url))  # Python 3
data=page.read().decode("utf8").split("<td><A HREF=")  # Python 3
tag=".fits"
endtag="</tr>"
for item in data:
    
    if ".fits" in item:
        start_pt = item.find("\"")
        end_pt = item.find("\"", start_pt + 1)
        url_want = item[start_pt + 1: end_pt]
        if url_want.find("spikes") == -1:
        
            if url_want[:4] == 'http':      
                filename0 = url_want
                arr_need = np.append(arr_need,filename0)
                fn_end = item.find(".fits\"")    
                
                if url_want.find("hmi.M") != -1:
                    fn_start = item.find("hmi.M")
                    filename = filename0[fn_start-1:fn_end]
                    
                    k = 0  # for adjusting
                    fyear = filename[10+k:14+k]
                    fmonth = filename[14+k:16+k]
                    fday = filename[16+k:18+k]
                    fhour = filename[19+k:21+k]
                    fmin = filename[21+k:23+k]
                    fsec = filename[23+k:25+k]
                    wavelength = 'magnetogram'
                    new_name = 'hmi_lev1_%s_%s_%s_%st%s_%s_%s_45z_image_lev1.fits' % (wavelength, fyear, fmonth, fday, fhour, fmin, fsec)  
                    arr_rename = np.append(arr_rename, new_name)
    
                elif url_want.find("hmi.I") != -1:
                    fn_start = item.find("hmi.I")
                    filename = filename0[fn_start-1:fn_end]
                    
                    k = 0  # for adjusting
                    fyear = filename[11+k:15+k]
                    fmonth = filename[15+k:17+k]
                    fday = filename[17+k:19+k]
                    fhour = filename[20+k:22+k]
                    fmin = filename[22+k:24+k]
                    fsec = filename[24+k:26+k]
                    wavelength = 'continuum'
                    new_name = 'hmi_lev1_%s_%s_%s_%st%s_%s_%s_45z_image_lev1.fits' % (wavelength, fyear, fmonth, fday, fhour, fmin, fsec) 
                    arr_rename = np.append(arr_rename, new_name)
                    
                elif url_want.find("aia.l") != -1:
                    fn_start = item.find("aia.l")
                    filename = filename0[fn_start-1:fn_end]
                    
                    if url_want.find("euv") != -1:
                        k = 1
                    else:
                        k = 0
                    
                    #k = 0  # for adjusting
                    fyear = filename[16+k:20+k]
                    fmonth = filename[21+k:23+k]
                    fday = filename[24+k:26+k]
                    fhour = filename[27+k:29+k]
                    fmin = filename[29+k:31+k]
                    fsec = filename[31+k:33+k]
                    wavelength_temp = filename[35+k:39+k]
                    if wavelength_temp[3] == '.':
                        wavelength = wavelength_temp[0:3]
                        cadence = 12
                    else:
                        wavelength = wavelength_temp
                        cadence = 24
                    
                    new_name = 'aia_lev1_%sa_%s_%s_%st%s_%s_%s_%iz_image_lev1.fits' % (wavelength, fyear, fmonth, fday, fhour, fmin, fsec, cadence) 
                    arr_rename = np.append(arr_rename, new_name)

#arr_need0 = arr_need[:8]  # for testing small number of files
#arr_rename0 = arr_rename[:8]  # for testing small number of files

arr_need0 = arr_need
arr_rename0 = arr_rename     

chunks = np.array_split(arr_need0, size)  # Split the data based on no. of processors
chunksB = np.array_split(arr_rename0, size)  # Split the data based on no. of processors
sub_arr = comm.scatter(chunks, root=0)
sub_arrB = comm.scatter(chunksB, root=0)
filesPart = get_data_fill(sub_arr, sub_arrB, directory)  # Do something with the array
filesTot = comm.gather(filesPart, root=0)  # Gather all the results  **check if this is ever different from actual, maybe just increasing counter even when miss file?

if rank == 0: 
    print("Downloaded %i/%i files" % (np.sum(filesTot), len(arr_rename0)))
#"""