# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 19:23:08 2017

@author: Brendan
"""

# maybe use this just for filling in.  time range seems to work just as well for initial pass - though maybe up to 5 min?  

from sunpy.net import vso
import astropy.units as u
import glob
import numpy as np

t1 = '2015/11/01 00:00:00'
t2 = '2015/11/01 00:01:00'
wavelength = 193
path_name = 'F:/Users/Brendan/Desktop/SolarProject/data2/20141115/193'
#path_name = '/mnt/data-solar/Gallagher/193_vso_test'


client=vso.VSOClient()  # establish connection to VSO database

qr=client.query(vso.attrs.Time(t1,t2), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
#res=client.get(qr, path='%s/{file}.fits' % path_name)  # leave the "{file}.fits" part alone  <-- ? 

#"""
arr_need = [] 
arr_need2 = [] 

for i in range(len(qr)):
    sec = qr[i].time[1][12:14]
    minute = qr[i].time[1][10:12]
    hour = qr[i].time[1][8:10]
    day = qr[i].time[1][6:8]
    month = qr[i].time[1][4:6]
    year = qr[i].time[1][0:4]
    time = '%s/%s/%s' ' %s:%s:%s' % (year, month, day, hour, minute, sec)
    #time2 = '%s/%s/%s' ' %s:%s:%s' % (year, month, day, hour, minute, sec+1)
    time = time.encode('utf8')
    #time2 = time2.encode('utf8')
    arr_need = np.append(arr_need, time)
    #arr_need2 = np.append(arr_need2, time2)
    
print qr
print arr_need

for i in range(len(arr_need)):
    
    #res=client.get(qr, path='%s/{file}.fits' % path_name).wait()  # leave the "{file}.fits" part alone  <-- ? 
    print arr_need[i]
    qr=client.query(vso.attrs.Time(arr_need[i],arr_need[i]), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
    print qr
    res=client.get(qr, path='%s/{file}.fits' % path_name).wait()  # leave the "{file}.fits" part alone  <-- ? 
    #print res
#"""