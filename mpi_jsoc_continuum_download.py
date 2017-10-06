# -*- coding: utf-8 -*-
"""
Created on Mon Oct 02 08:55:48 2017

@author: Brendan
"""

"""
###########################
## Downloads .FITS files from JSOC temporary url
###########################
"""

import numpy as np
from mpi4py import MPI
import urllib
import urllib2

#"""
def get_data_fill(arr_need, arr_rename, directory):
    
    if rank == 0:
        print "Please wait while request is being processed."    
        sub_arr_len = len(arr_need)
    
    counter = 0
        
    # loop through the array of needed files, requesting them one at a time         
    for i in range(len(arr_need)):
    
        urllib.urlretrieve("%s" % arr_need[i], "%s/FITS/%s/%s/%s" % (directory, date, 'hmi_continuum', arr_rename[i]))
        
        counter += 1
        
        if rank == 0:
            print "Probably downloading file %i/%i" % (counter*size, sub_arr_len*size)
    return counter


              
comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"

size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

import sys

jsoc_url = sys.argv[1]
directory = sys.argv[2]
date = sys.argv[3]


arr_need = []
arr_rename = []

page=urllib2.urlopen(jsoc_url)
data=page.read().split("<td><A HREF=")
tag=".fits"
endtag="</tr>"
for item in data:
    if ".fits" in item:
        
        if item[1:5] == 'http':
                
                filename0 = item[1:97]
                arr_need = np.append(arr_need,filename0)
                filename = filename0[49:]
                
                k = -3  # for adjusting
                fyear = filename[14+k:18+k]
                fmonth = filename[18+k:20+k]
                fday = filename[20+k:22+k]
                fhour = filename[23+k:25+k]
                fmin = filename[25+k:27+k]
                fsec = filename[27+k:29+k]
                #wavelength = filename[33+k:37+k]  # maybe use - characters back from end - then "if 1600/1700 do this__"
                new_name = 'hmi_lev1_continuum_%s_%s_%st%s_%s_%s_45z_image_lev1.fits' % (fyear, fmonth, fday, fhour, fmin, fsec) 
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
    print "Downloaded %i/%i files" % (np.sum(filesTot), len(arr_rename0))
#"""