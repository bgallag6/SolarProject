# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 22:33:10 2016

@author: Brendan
"""
### Example full-module function calls 


import numpy as np
import SolSpec as ss

"""
** not in current module **
## download data
"""
#r = ss.get_data(wavelength=1600, time_begin='2013/06/26 00:00:00', time_end='2013/06/26 00:05:00', path_name='F:/Users/Brendan/Desktop/SolarProject/data2/20130626/1600')


"""
** not in current module **
## download data (fill in missing)
"""
#r = ss.get_data_fill(wavelength=193, time_begin='2013/05/30 00:00:00', time_end='2013/05/30 11:59:59', path_name='F:/Users/Brendan/Desktop/SolarProject/FITS/20130530/193')


"""
## arcsecond to pixel + subregion
"""
#x1 = -550
#x2 = 550
#y1 = -500
#y2 = 500

#x1 = 1850
#x2 = 3050
#y1 = 2000
#y2 = 2900

#ss.arc2pix(x1,x2,y1,y2, image = 'F:/Users/Brendan/Desktop/SolarProject/FITS/20140818/1600/aia_lev1_1600a_2014_08_18t18_00_16_12z_image_lev1.fits')
#ss.arc2pix(x1,x2,y1,y2, image = 'S:/FITS/20130626/193/aia_lev1_193a_2013_06_26t06_00_06_84z_image_lev1.fits')

#"""
#directory = 'F:/Users/Brendan/Desktop/SolarProject'
#directory = '/mnt/data/Gallagher'
#date = '20160426'
#wavelength = 1600
#"""
import sys
directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])

"""
## create derotated region datacube
"""
#ss.datacube(directory='F:/Users/Brendan/Desktop/SolarProject/data/20120923/304', date='20120923', wavelength=304, sub_reg_coords=[-528,-132,-100,100], coords_type='arc', bin_frac=1)
ss.datacube(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength, sub_reg_coords=[-395,395,-375,375], coords_type='arc', bin_frac=1)


"""
## fft-averaging + 3x3 averaging
"""
#DATA = np.load('%s/DATA/Temp/%s/%i/*rebin1.npy')
#TIME = np.load('%s/DATA/Temp/%s/%i/*time.npy')
#EXPOSURE = np.load('%s/DATA/Temp/%s/%i/*exposure.npy')

ss.fft_avg(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength, num_seg = 6)  # 3 seg for 20120923, 6 seg for rest
#spectra_array = ss.fft_avg(directory='%s' % (directory), date='%s' % (date), wavelength= '%i' % (wavelength), datacube = DATA, timeseries = TIME, exposure_array = EXPOSURE, num_seg = 6)
#np.save('F:/Users/Brendan/Desktop/SolarProject/data/20120923/304/20120923_304_-100_100i_-528_-132j_spectra', spectra_array)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_193_2300_2600_2200_3000', spectra)  # now this


"""
## fft-averaging + 3x3 averaging (overlapping time-segments)
"""
#ss.fft_overlap(directory='%s' % directory, date='%s' % date, wavelength=wavelength, window_length='02:00', overlap_pct=75, pixel_box=True)


"""
## memory mapping 
"""

ss.mem_map(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength)
#spectra_array = ss.fft_avg(directory='%s' % (directory), date='%s' % (date), wavelength= '%i' % (wavelength), datacube = DATA, timeseries = TIME, exposure_array = EXPOSURE, num_seg = 6)
#np.save('F:/Users/Brendan/Desktop/SolarProject/data/20120923/304/20120923_304_-100_100i_-528_-132j_spectra', spectra_array)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_193_2300_2600_2200_3000', spectra)  # now this
"""


"""
## generate heatmaps
"""
#ss.heatmap(directory= '%s' % (directory), date='%s' % (date), wavelength= '%i' % (wavelength))
"""