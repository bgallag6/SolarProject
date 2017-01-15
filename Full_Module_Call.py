# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 22:33:10 2016

@author: Brendan
"""
### Example full-module function calls 

import numpy as np
import Full_Module_Automated as fm
import h5py

"""
## download data
"""
#r = fm.get_data(wavelength=193, time_begin='2014/09/02 00:00:00', time_end='2014/09/02 12:00:00', path_name='F:/Users/Brendan/Desktop/SolarProject/data/20140902/193')


"""
## download data (fill in missing)
"""
#r = fm.get_data_fill(wavelength=193, cadence=12, time_begin='2014/09/02 00:00:00', time_end='2014/09/02 12:00:00', path_name='F:/Users/Brendan/Desktop/SolarProject/data/20140902/193')


"""
## generate heatmaps
"""
#with h5py.File("D:/Users/Brendan/Desktop/SolarProject/hdf5/20141025_1600_1350_1650i_2150_3150j_float_dogbox_24hr.hdf5",'r') as f:
    
    #HEATMAPS = np.array(f['params'])
    #VISUAL = np.array(f['visual'])
    #r = fm.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20141025', wavelength=1600, path_name='C:/Users/Brendan/Desktop/PDFs')
    
#HEATMAPS = np.load('C:/Users/Brendan/Desktop/SDO/param_20120923_304A_(528)_(132)x_(100)_100y.npy')
#HEATMAPS = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20120923_193A_(528)_(132)x_(100)_100y_numpy_dogbox.npy')
#VISUAL = np.load('F:/Users/Brendan/Desktop/SolarProject/visual/visual_20120923_304A_(528)_(132)x_(100)_100y_float.npy')
#r = fm.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20120923', wavelength=304, path_name='C:/Users/Brendan/Desktop/PHYS 326')



"""
## arcsecond to pixel + subregion
"""
#x1 = -100
#x2 = 600
#y1 = 0
#y2 = 500

#x1 = 1850
#x2 = 3050
#y1 = 2000
#y2 = 2900

#fm.arc2pix(x1,x2,y1,y2, image = 'F:/Users/Brendan/Desktop/SolarProject/data/20140902/193/aia_lev1_193a_2014_09_02t05_59_54_84z_image_lev1.fits.fits')
#fm.pix2arc(x1,x2,y1,y2, image = 'F:/Users/Brendan/Desktop/SolarProject/data/20140902/193/aia_lev1_193a_2014_09_02t05_59_54_84z_image_lev1.fits.fits')


"""
## create derotated region datacube
"""
#fm.datacube(directory='F:/Users/Brendan/Desktop/SolarProject/data/20140902/193', date='20140902', wavelength=193, sub_reg_coords=[1850,3050,2000,2900], coords_type='pix', bin_frac=2)



"""
## fft-averaging
"""
#DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/datacubes/SDO_20120923_304A_(528)_(132)x_(100)_100y_float.npy')
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/time_arrays/SDO_20120923_304A_(528)_(132)x_(100)_100y_float_time.npy')

#spectra_array = fm.fft_avg(datacube = DATA, timeseries = TIME, num_seg = 3)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_array_FFT3_20120923_304A_(528)_(132)x_(100)_100y.npy', spectra_array)


"""
## 3x3 pixel box averaging + fitting
"""
#SPECTRA = np.load('C:/Users/Brendan/Desktop/SDO/spectra_array_FFT6_20130530_1600_2300_2600i_2200_3000j_rebin4_rev_t_interp.npy')
SPECTRA = np.load('F:/Users/Brendan/Desktop/SolarProject/spectra_array_FFT/spectra_array_FFT6_20130530_193_2300_2600_2200_3000_float_numpy.npy')
#SPECTRA = spectra_array
params, spectra, M2_fit = fm.spec_fit(spectra_array = SPECTRA)

np.save('C:/Users/Brendan/Desktop/SDO/param_20130530_193_2300_2600_2200_3000', params)
np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_193_2300_2600_2200_3000', spectra)
np.save('C:/Users/Brendan/Desktop/SDO/M2_20130530_193_2300_2600_2200_3000', M2_fit)
#np.save('C:/Users/Brendan/Desktop/SDO/uncertainties_20130815_193_1000_1600i_1950_2950j_rebin2', Uncertainties)  # if want to keep?
#np.save('C:/Users/Brendan/Desktop/SDO/diffm1m2_20130815_193_1000_1600i_1950_2950j_rebin2', diffM1M2)  # if want to keep?


"""
## parameter masking  (doesn't work - wont stay connected to plot as function call)
"""
#HEATMAPS = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20130530_1600_0_296_0_634_numpy.npy')
#fm.mask_param(heatmaps = HEATMAPS)