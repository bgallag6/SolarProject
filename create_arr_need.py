# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 05:02:58 2017

@author: Brendan
"""

# loop through flist couple times until get all.

from sunpy.net import vso
import astropy.units as u
import glob
import numpy as np

import sys

directory = sys.argv[1]
date_slash = sys.argv[2]
date_reg = sys.argv[3]
time_begin = sys.argv[4]
time_end = sys.argv[5]
wavelength = int(sys.argv[6])

print "Please wait while request is being processed."

client=vso.VSOClient()  # establish connection to VSO database

# extract yyyy/mm/dd, hour, minute, and second values from start of time-range
Y1 = str(date_slash)   
H1 = int(time_begin[0:2])
M1 = int(time_begin[3:5])
S1 = int(time_begin[6:8])

# extract yyyy/mm/dd, hour, minute, and second values from end of time-range
Y2 = str(date_slash)   
H2 = int(time_end[0:2])
M2 = int(time_end[3:5])
S2 = int(time_end[6:8])

# create time strings to pass to initial query request
T1 = ('%-11s''%02d'':''%02d'':''%02d' % (Y1,H1,M1,S1))
T2 = ('%-11s''%02d'':''%02d'':''%02d' % (Y2,H2,M2,S2))

# query request to determine total number of files in time-range
qr=client.query(vso.attrs.Time(T1,T2), vso.attrs.Instrument('aia'), vso.attrs.Wavelength(wavelength * u.AA, wavelength * u.AA))

arr_all = [] 

for i in range(len(qr)):
    sec = qr[i].time[1][12:14]
    minute = qr[i].time[1][10:12]
    hour = qr[i].time[1][8:10]
    day = qr[i].time[1][6:8]
    month = qr[i].time[1][4:6]
    year = qr[i].time[1][0:4]
    time = '%s/%s/%s' ' %s:%s:%s' % (year, month, day, hour, minute, sec)
    time = time.encode('utf8')
    arr_all = np.append(arr_all, time)

num_files = len(qr)
print num_files
#cadence = cadence  # set cadence to specified value *took out because extracting times from query
    
flist = glob.glob('%s/FITS/%s/%i/*.fits' % (directory, date_reg, wavelength))

l = len(flist)
   
# create searchable array of images that have already been downloaded
    
adj = 0  # adjust 5 characters for 20130815 193 dataset (doesn't have extra '.fits')
#adj = -5  # for the datasets left that have extra '.fits.

arr_have = []

if l > 0: 
    l_fname = len(flist[0])
    for i in range(0,l):
        x = flist[i]
        h = int(x[(l_fname-28+adj):(l_fname-26+adj)])
        m = int(x[(l_fname-25+adj):(l_fname-23+adj)])
        s = int(x[(l_fname-22+adj):(l_fname-20+adj)])        
        t = ('%-11s''%02d'':''%02d'':''%02d' % (Y1,h,m,s))
        arr_have.append(t)
#print arr_have
print len(arr_have)    
#print arr_have[0]
#print arr_all[0]

# compare array_all to array_have to determine array_need
arr_need = []   
for i in range(0,num_files):  # changed from range(0,l)
    z = arr_all[i] in arr_have    
    if z == False:
        arr_need.append(arr_all[i])
#print arr_need
#print len(arr_need)    

print ""
print "Query request completed, %d files are needed." % len(arr_need)

np.save('%s/FITS/%s/%i/arr_need.npy' % (directory, date_reg, wavelength), arr_need)