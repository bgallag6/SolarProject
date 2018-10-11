# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 17:26:53 2018

@author: Brendan
"""

import numpy as np
import glob
import astropy.units as u
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.time import Time

"""
# load .fits, select subregion, save trimmed image

in_dir = 'S:/FITS/20120606/1600_middle/'

flist = sorted(glob.glob('%s*' % in_dir))

aia0 = sunpy.map.Map(flist[0])

bl = SkyCoord(125*u.arcsec, 200*u.arcsec, frame=aia0.coordinate_frame)  
tr = SkyCoord(300*u.arcsec, 300*u.arcsec, frame=aia0.coordinate_frame) 

out_dir = 'C:/Users/Brendan/Desktop/test_data'

for filename in flist:
    aia = sunpy.map.Map(filename)
    sub_map = aia.submap(bl, tr)
    sub_map.save('%s/%s' % (out_dir, filename[len(in_dir):]))
"""


"""
# verify the timestamps of the new files

direct = 'C:/Users/Brendan/Desktop/test_data/*'

flist = sorted(glob.glob('%s' % direct))
nf = len(flist)

timestamp = np.zeros((nf))  # seconds elapsed since first image

counter = 0

for filename in flist:
    aia = sunpy.map.Map(filename)
    timestamp[counter] = Time(aia.date).jd  # extract julian day time from each image
    counter += 1
    
timestamp -= timestamp[0]  # calculate time since first image
timestamp = np.around(timestamp*86400)  # get the time value in seconds, and round to nearest whole number
"""





"""
#y = aia.submap(bl, tr).data   
y = aia.submap(bl, tr)
#e = aia.exposure_time
#t = aia.date
#print(e, t)

out_dir = 'C:/Users/Brendan/Desktop/test_data'

y.save('%s/%s' % (out_dir, fname))
"""


#sunpy.map.Map.save("C:/Users/Brendan/Desktop/sunpy_try1.fits")
"""
aia2 = sunpy.map.Map('%s/%s' % (out_dir, fname))
y2 = aia2.data
e2 = aia2.exposure_time
t2 = aia2.date
print(e2, t2)
"""