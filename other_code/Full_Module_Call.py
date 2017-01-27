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
#r = fm.get_data_fill(wavelength=171, cadence=12, time_begin='2013/05/30 00:00:00', time_end='2013/05/30 12:00:00', path_name='F:/Users/Brendan/Desktop/SolarProject/data/20130530/171')


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
#fm.datacube(directory='F:/Users/Brendan/Desktop/SolarProject/data/20130530/1700', date='20130530', wavelength=1700, sub_reg_coords=[2200,3000,2300,2600], coords_type='pix', bin_frac=1)



"""
## fft-averaging + 3x3 averaging
"""
#DATA = np.load('F:/Users/Brendan/Desktop/SolarProject/datacubes/20130530_1700_2300_2600i_2200_3000j_data_rebin1.npy')
#TIME = np.load('F:/Users/Brendan/Desktop/SolarProject/time_arrays/20130530_1700_2300_2600i_2200_3000j_time.npy')

#spectra_array = fm.fft_avg(datacube = DATA, timeseries = TIME, num_seg = 6)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_1700_2300_2600i_2200_3000j', spectra_array)
#np.save('C:/Users/Brendan/Desktop/SDO/spectra_20130530_193_2300_2600_2200_3000', spectra)  # now this


"""
## spectra fitting
"""
#SPECTRA = np.load('C:/Users/Brendan/Desktop/SDO/spectra_20130530_1600_2300_2600i_2200_3000j_data_rebin4.npy')
#SPECTRA = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/spectra_20130815_193_1000_1600i_1950_2950j_rebin2.npy')
#SPECTRA = spectra_array
#params, M2_fit = fm.spec_fit(spectra_array = SPECTRA)

#np.save('C:/Users/Brendan/Desktop/SDO/param_20130530_1600_2300_2600i_2200_3000j_data_rebin4_float_32', params)
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
    
#HEATMAPS = np.load('C:/Users/Brendan/Desktop/SDO/param_20130530_1600_2300_2600i_2200_3000j_data_rebin4_float_32.npy')
#HEATMAPS = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20130530_193_-500_500i_-500_500j.npy')
#VISUAL = np.load('F:/Users/Brendan/Desktop/SolarProject/visual/visual_20130530_193_-500_500i_-500_500j.npy')
#r = fm.heatmap(heatmaps = HEATMAPS, visual = VISUAL, date = '20130530', wavelength=193, path_name='C:/Users/Brendan/Desktop/PHYS 326/full_region_aspect2')


"""
## parameter masking  (doesn't work - wont stay connected to plot as function call)
"""
#HEATMAPS = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20130530_1600_0_296_0_634_numpy.npy')
#fm.mask_param(heatmaps = HEATMAPS)