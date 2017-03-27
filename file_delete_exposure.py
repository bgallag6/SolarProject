# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 14:19:53 2017

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
import numpy as np
import astropy.units as u
import os

#bad = Map('F:/Users/Brendan/Desktop/SolarProject/FITS/20140317/304/aia_lev1_304a_2014_03_17t21_00_23_51z_image_lev1.fits')
#good = Map('F:/Users/Brendan/Desktop/SolarProject/FITS/20140317/304/aia_lev1_304a_2014_03_17t21_00_43_12z_image_lev1.fits')

#print bad.exposure_time
#print good.exposure_time


#plt.figure()
#plt.imshow(bad.data)

#plt.figure()
#plt.imshow(good.data)

flist = glob.glob('F:/Users/Brendan/Desktop/SolarProject/FITS/20160905/1712/aia*.fits')

nf = len(flist)

ex = []


for filename in flist:
    ex.append(Map(filename).exposure_time.value)
    
med = np.median(ex)
med_top = med*1.5
med_bot = med*0.5

flist_del = []

count = 0
num_del = 0
for filename in flist:
    if ex[count] > med_top or ex[count] < med_bot:
        flist_del.append(filename)
        num_del += 1
    count += 1
       
print("%i files will be deleted." % num_del)
raw_input("Press [Enter] to continue...")

for filename in flist_del:
        os.remove("%s" % filename)
        print("Deleted file: %s" % filename)