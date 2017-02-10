# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 22:33:10 2016

@author: Brendan
"""
### Example full-module function calls 

import numpy as np
import SolSpec as ss
import h5py

"""
## download data
"""
#r = ss.get_data(wavelength=304, time_begin='2014/11/15 00:00:00', time_end='2014/11/15 12:00:00', path_name='F:/Users/Brendan/Desktop/SolarProject/data/20141115/304')


"""
## download data (fill in missing)
"""
#r = ss.get_data_fill(wavelength=304, cadence=12, time_begin='2014/11/15 01:40:00', time_end='2014/11/15 12:00:00', path_name='F:/Users/Brendan/Desktop/SolarProject/data/20141115/304')


"""
## arcsecond to pixel + subregion
"""
#x1 = -100
#x2 = 100
#y1 = -100
#y2 = 100

#x1 = 1850
#x2 = 3050
#y1 = 2000
#y2 = 2900

#ss.arc2pix(x1,x2,y1,y2, image = 'F:/Users/Brendan/Desktop/SolarProject/data/20130626/304/aia_lev1_304a_2013_06_26t00_00_07_12z_image_lev1.fits.fits')
#ss.pix2arc(x1,x2,y1,y2, image = 'F:/Users/Brendan/Desktop/SolarProject/data/20140902/193/aia_lev1_193a_2014_09_02t05_59_54_84z_image_lev1.fits.fits')


"""
## create derotated region datacube (old)
"""
#ss.datacube(directory='F:/Users/Brendan/Desktop/SolarProject/data/20130530/193', date='20130530', wavelength=193, sub_reg_coords=[2200,3000,2300,2600], coords_type='pix', bin_frac=1)


"""
## create derotated region datacube (int - new)
"""
#ss.datacube_int(directory='F:/Users/Brendan/Desktop/SolarProject/data/20130626/335', date='20130626', wavelength=335, sub_reg_coords=[-180,0,170,375], coords_type='arc', bin_frac=1)


"""
## fft-averaging + 3x3 averaging (old)
"""
#DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/datacubes/20130530_1700_2300_2600i_2200_3000j_data_rebin1.npy')
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/time_arrays/20130530_1700_2300_2600i_2200_3000j_time.npy')

#spectra_array = ss.fft_avg(datacube = DATA, timeseries = TIME, num_seg = 6)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_1700_2300_2600i_2200_3000j', spectra_array)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_193_2300_2600_2200_3000', spectra)  # now this


"""
## fft-averaging + 3x3 averaging (int - new)
"""
#DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/335/20130626_335_170_375i_-180_0j_data_rebin1.npy')
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/335/20130626_335_170_375i_-180_0j_time.npy')
#EXPOSURE = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/335/20130626_335_170_375i_-180_0j_exposure.npy')

#spectra_array = ss.fft_avg_int(datacube = DATA, timeseries = TIME, exposure_array = EXPOSURE, num_seg = 6)
#np.save('F:/Users/Brendan/Desktop/SolarProject/data/20130626/335/20130626_335_170_375i_-180_0j_spectra', spectra_array)
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
    
#HEATMAPS = np.load('C:/Users/Brendan/Desktop/SDO/20130530_193_2300_2600i_2200_3000j_rebin1_params_mpi.npy')
#HEATMAPS = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/304/20130626_304_-500_500i_-500_500j_param.npy')
#VISUAL = np.load('F:/Users/Brendan/Desktop/SolarProject/data/20130626/304/20130626_304_-500_500i_-500_500j_visual.npy')
#r = ss.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20130626', wavelength=304, path_name='C:/Users/Brendan/Desktop/20130626/')
