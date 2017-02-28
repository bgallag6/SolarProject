# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 22:33:10 2016

@author: Brendan
"""
### Example full-module function calls 

import numpy as np
import SolSpec as ss
#import h5py

"""
## download data
"""
#r = ss.get_data(wavelength=1600, time_begin='2013/06/26 08:23:00', time_end='2013/06/26 12:00:00', path_name='F:/Users/Brendan/Desktop/SolarProject/data/20130626/1600')


"""
## download data (fill in missing)
"""
#r = ss.get_data_fill(wavelength=193, cadence=12, time_begin='2012/09/23 00:00:00', time_end='2012/09/23 06:00:00', path_name='F:/Users/Brendan/Desktop/SolarProject/data/20120923/193')


"""
## arcsecond to pixel + subregion
"""
#x1 = -450
#x2 = 550
#y1 = -500
#y2 = 500

#x1 = 1850
#x2 = 3050
#y1 = 2000
#y2 = 2900

#ss.arc2pix(x1,x2,y1,y2, image = 'F:/Users/Brendan/Desktop/SolarProject/data/20130626/171/aia_lev1_171a_2013_06_26t00_00_11_34z_image_lev1.fits.fits')
#ss.pix2arc(x1,x2,y1,y2, image = 'F:/Users/Brendan/Desktop/SolarProject/data/20140902/193/aia_lev1_193a_2014_09_02t05_59_54_84z_image_lev1.fits.fits')

#directory = '/mnt/data/Gallagher'
directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20120923'
wavelength = 94

"""
## create derotated region datacube
"""
#ss.datacube(directory='F:/Users/Brendan/Desktop/SolarProject/data/20120923/304', date='20120923', wavelength=304, sub_reg_coords=[-528,-132,-100,100], coords_type='arc', bin_frac=1)
ss.datacube(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength, sub_reg_coords=[-528,-132,-100,100], coords_type='arc', bin_frac=1)
#ss.datacube(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength, sub_reg_coords=[2200,3000,2300,2600], coords_type='pix', bin_frac=1)


"""
## fft-averaging + 3x3 averaging
"""
#DATA = np.load('%s/DATA/Temp/%s/%i/*rebin1.npy')
#TIME = np.load('%s/DATA/Temp/%s/%i/*time.npy')
#EXPOSURE = np.load('%s/DATA/Temp/%s/%i/*exposure.npy')

ss.fft_avg(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength, num_seg = 3)
#spectra_array = ss.fft_avg(directory='%s' % (directory), date='%s' % (date), wavelength= '%i' % (wavelength), datacube = DATA, timeseries = TIME, exposure_array = EXPOSURE, num_seg = 6)
#np.save('F:/Users/Brendan/Desktop/SolarProject/data/20120923/304/20120923_304_-100_100i_-528_-132j_spectra', spectra_array)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_193_2300_2600_2200_3000', spectra)  # now this


"""
## memory mapping 
"""
#DATA = np.load('%s/DATA/Temp/%s/%i/*rebin1.npy')
#TIME = np.load('%s/DATA/Temp/%s/%i/*time.npy')
#EXPOSURE = np.load('%s/DATA/Temp/%s/%i/*exposure.npy')

ss.mem_map(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength)
#spectra_array = ss.fft_avg(directory='%s' % (directory), date='%s' % (date), wavelength= '%i' % (wavelength), datacube = DATA, timeseries = TIME, exposure_array = EXPOSURE, num_seg = 6)
#np.save('F:/Users/Brendan/Desktop/SolarProject/data/20120923/304/20120923_304_-100_100i_-528_-132j_spectra', spectra_array)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_193_2300_2600_2200_3000', spectra)  # now this


"""
## spectra fitting
"""
#SPECTRA = np.load('C:/Users/Brendan/Desktop/SDO/20130530_193_2300_2600i_2200_3000j_rebin1_spectra_mpi.npy')
#SPECTRA = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/335/20130626_335_170_375i_-180_0j_spectra.npy')
#SPECTRA = spectra_array
#params, M2_fit = ss.spec_fit(spectra_array = SPECTRA)

#np.save('F:/Users/Brendan/Desktop/SolarProject/data/20130626/335/20130626_335_170_375i_-180_0j_param', params)
#np.save('C:/Users/Brendan/Desktop/SDO/M2_20130530_1700_2300_2600i_2200_3000j', M2_fit)
#np.save('C:/Users/Brendan/Desktop/SDO/uncertainties_20130815_193_1000_1600i_1950_2950j_rebin2', Uncertainties)  # if want to keep?
#np.save('C:/Users/Brendan/Desktop/SDO/diffm1m2_20130815_193_1000_1600i_1950_2950j_rebin2', diffM1M2)  # if want to keep?


"""
## generate heatmaps
"""
#with h5py.File("D:/Users/Brendan/Desktop/SolarProject/hdf5/20141025_1600_1350_1650i_2150_3150j_float_dogbox_24hr.hdf5",'r') as f:
    
    #HEATMAPS = np.array(f['params'])
    #VISUAL = np.array(f['visual'])
    #r = fm.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20141025', wavelength=1600, path_name='C:/Users/Brendan/Desktop/PDFs')
    
#HEATMAPS = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_param.npy')
#HEATMAPS = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20120923/171/20120923_171_-100_100i_-528_-132j_param.npy')
#VISUAL = np.load('C:/Users/Brendan/Desktop/20130626_193_-500_500i_-500_600j_visual.npy')
#VISUAL = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20120923/171/20120923_171_-100_100i_-528_-132j_visual.npy')
#ss.heatmap(directory= '%s' % (directory), date='%s' % (date), wavelength= '%i' % (wavelength))
#r = ss.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20130626', wavelength=193, path_name='C:/Users/Brendan/Desktop/test_delete/')
