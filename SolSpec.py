# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 22:30:43 2016
@author: Brendan Gallagher

This module can be used to:
1. Download .fits image files
2. Extract subregion and deroate datacube
3. Calculate the time-segment-averaged FFT for each pixel in region
4. Calculate a 3x3 pixel-box average for each pixel
5. Fit the resulting FFT spectra to models M1 and M2
6. Generate 'heatmaps' for each of the parameters of M2

"""

# maybe make separate function with whole program -- from get_data to heatmaps?

# maybe put module imports inside function definitions (so not executing all for each)
# delete unused modules for each function

# maybe write script to compare results easier

# update 1/21:
# maybe make duplicate fft_avg function, one for computers that have accelerate, one for don't

# as of 1/23: have tested all functions to be working (except parameter masking tool)

"""
############################
############################
# download data 
############################
############################
"""

# where should I specify path to save (within function or before/during/after function call)
# put in error handling if start date is after end date? also if data is entered incorrectly - not yyyy/mm/dd format

# generalize for 1600 - string bounds wont work - add 1

# maybe make first past be much quicker - like 10 min at a time

# for 20130815 193 - files are just ".fits" not ".fits.fits" - screws up call obviously

from sunpy.net import vso
import astropy.units as u
import glob
import numpy as np

def get_data(wavelength, time_begin, time_end, path_name):
    """
    Downloads .fits image files from database. 
       
    wavelength : 
        Wavelength to be downloaded. (Integer)
       
    time_begin : 
        Beginning of time range in YYYY/MM/DD HH:MM:SS format. (String)
     
    time_end : 
        Ending of time range in YYYY/MM/DD HH:MM:SS format. (String)
                   
    path_name : 
        The directory to which the files should be downloaded. (String)
      
    Example:
    ::
        ss.get_data(wavelength=1600, time_begin='2016/09/23 00:00:00', time_end='2016/09/23 00:05:00',
           path_name='C:/Users/Brendan/Desktop/SDO_test')
    """
    
    print ""
    print "Please wait while request is being processed."

    client=vso.VSOClient()  # establish connection to VSO database
    
    # extract yyyy/mm/dd, hour, minute, and second values from start of time-range
    Y1 = str(time_begin[0:11])   
    H1 = int(time_begin[11:13])
    M1 = int(time_begin[14:16])
    S1 = int(time_begin[17:19])
    
    # extract yyyy/mm/dd, hour, minute, and second values from end of time-range
    Y2 = str(time_end[0:11])   
    H2 = int(time_end[11:13])
    M2 = int(time_end[14:16])
    S2 = int(time_end[17:19])
    
    # create time strings to pass to initial query request
    T1 = ('%s''%02d'':''%02d'':''%02d' % (Y1,H1,M1,S1))
    T2 = ('%s''%02d'':''%02d'':''%02d' % (Y2,H2,M2,S2))
    
    # query request to determine total number of files in time-range
    qr=client.query(vso.attrs.Time(T1,T2), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
    
    num_files = len(qr)  # total number of files in time-range
    duration = (H2-H1)*3600 + (M2-M1)*60  # total length of time-range
    spacing = duration / num_files  # determine average image cadence
    if spacing > 10 and spacing < 14:
        cadence = 12
    elif spacing > 22 and spacing < 26:
        cadence = 24
    
    print ""    
    print "This request will download the dataset of wavelength %d" % wavelength
    print "in the timerange from %s to %s." % (T1,T2)
    print "Image files will be saved in the directory: %s" % path_name
    print ""
    print "This will download %d files at %d-second cadence." % (num_files,cadence)
    print ""
    raw_input("Press [Enter] to continue...")
    
    # set starting values of iteration to specified start time
    s1 = S1
    s2 = s1
    m1 = M1
    m2 = m1 + 1
    #m2 = m1 + mins  # if want to try for faster 
    h1 = H1
    h2 = h1
    
    # loop through entire time range, 1 minute at a time (so as not to overload the server with requests)
    for i in range(0,(duration/60)-1):
    #for i in range(0,(duration/(60*mins))-1):  # try for quicker download using custom minutes-per-query
        m1 += 1
        m2 += 1
        #m1 += mins
        #m2 += mins
        if m1 >= 60:
            m1 -= 60
            h1 += 1
        if m2 >= 60:
            m2 -= 60
            h2 += 1
        t1 = ('%s''%02d'':''%02d'':''%02d' % (Y1,h1,m1,s1))
        print t1
        t2 = ('%s''%02d'':''%02d'':''%02d' % (Y1,h2,m2,s2))
        print t2
        
        qr=client.query(vso.attrs.Time(t1,t2), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
        #print qr
        res=client.get(qr, path='%s/{file}.fits' % path_name).wait()  # leave the "{file}.fits" part alone  <-- ? 
        #print res
    
    ## after initial passthrough, determine files downloaded and those that still need to be    
        
    # get list of files downloaded    
    flist = glob.glob('%s/*.fits' % path_name)
    
    l = len(flist)
     
    l_fname = len(flist[0])
   
    # find first file after 00:00:00 - set time base to that
    
    # loop through flist couple times until get all.
    
    arr_have = []
    
    # create searchable array of images that have already been downloaded
    for i in range(0,l):
        x = flist[i]
        h = int(x[(l_fname-33):(l_fname-31)])
        m = int(x[(l_fname-30):(l_fname-28)])
        s = int(x[(l_fname-27):(l_fname-25)])
        t = ('%s''%02d'':''%02d'':''%02d' % (Y1,h,m,s))
        arr_have.append(t)
    #print arr_have
    
    # using the first image file as the starting time, add determined cadence to generate searchable array of all possible files
    # adding ' -1' to m3 because initial download starts at first minute
    f_first = flist[0]
    h3 = int(f_first[(l_fname-33):(l_fname-31)])
    m3 = int(f_first[(l_fname-30):(l_fname-28)]) - 1 
    s3 = (60 +  (int(f_first[(l_fname-27):(l_fname-25)])) - (int(np.floor(60/cadence))*cadence))
    #t_first = ('%s''%02d'':''%02d'':''%02d' % (Y1,H1,(m_first - 1),s_first))
       
    
    # create array of all possible files
    arr_all = []
    for i in range(0,num_files):  # changed from range(0,l)
        t3 = ('%s''%02d'':''%02d'':''%02d' % (Y1,h3,m3,s3))
        arr_all.append(t3)
        s3 += cadence
        if s3 >= 60:
            s3 -= 60
            m3 += 1
            if m3 >= 60:
                m3 -= 60
                h3 += 1 
    #print arr_all
    
    
    # compare array_all to array_have to determine array_need
    arr_need = []
    for i in range(0,num_files):  # changed from range(0,l)
        z = arr_all[i] in arr_have    
        if z == False:
            arr_need.append(arr_all[i])
    #print arr_need
    
    print ""
    print "After the initial pass, still need %d files." % len(arr_need)
    
    # loop through the array of needed files, requesting them one at a time        
    for i in range(0,len(arr_need)):
        qr=client.query(vso.attrs.Time(arr_need[i],arr_need[i]), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
        res=client.get(qr, path='%s/{file}.fits' % path_name).wait()
  
      

"""
############################
############################
# re-download data 
############################
############################
"""        
        
### maybe put as an if statement in the original - or specify argument in original call
# like (first / new   vs re / fill)
        
# prof weigel mentioned that could instead use push/pop to move first element of array to bottom
# so that program doesn't get stuck on one file, also that it would be better to use the ISO module for date/time
        
# when redownloading - takes first file and subtracts one minute - so have to delete first minute of files

        
from sunpy.net import vso
import astropy.units as u
import glob
import numpy as np

def get_data_fill(wavelength, cadence, time_begin, time_end, path_name):
    """
    Downloads .fits image files from database. 
    
    wavelength : 
        Wavelength to be downloaded. (Integer)
       
    time_begin : 
        Beginning of time range in YYYY/MM/DD HH:MM:SS format. (String)
     
    time_end : 
        Ending of time range in YYYY/MM/DD HH:MM:SS format. (String)
                   
    path_name : 
        The directory to which the files should be downloaded. (String)
      
    Example:
    ::
        ss.get_data_fill(wavelength=171, cadence=12, time_begin='2013/08/15 00:00:00',
           time_end='2013/08/15 12:00:00', path_name='F:/SDO/data/20130815/171')
    """
    
    print ""
    print "Please wait while request is being processed."

    client=vso.VSOClient()  # establish connection to VSO database
    
    # extract yyyy/mm/dd, hour, minute, and second values from start of time-range
    Y1 = str(time_begin[0:11])   
    H1 = int(time_begin[11:13])
    M1 = int(time_begin[14:16])
    S1 = int(time_begin[17:19])
    
    # extract yyyy/mm/dd, hour, minute, and second values from end of time-range
    Y2 = str(time_end[0:11])   
    H2 = int(time_end[11:13])
    M2 = int(time_end[14:16])
    S2 = int(time_end[17:19])
    
    # create time strings to pass to initial query request
    T1 = ('%s''%02d'':''%02d'':''%02d' % (Y1,H1,M1,S1))
    T2 = ('%s''%02d'':''%02d'':''%02d' % (Y2,H2,M2,S2))
    
    # query request to determine total number of files in time-range
    qr=client.query(vso.attrs.Time(T1,T2), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
    
    num_files = len(qr)
    
    cadence = cadence  # set cadence to specified value
        
    flist = glob.glob('%s/*.fits' % path_name)
    
    l = len(flist)
     
    l_fname = len(flist[0])
   
    # find first file after 00:00:00 - set time base to that
    
    # loop through flist couple times until get all.
    
    
    # create searchable array of images that have already been downloaded
        
    #adj = 5  # adjust 5 characters for 20130815 193 dataset (doesn't have extra '.fits')
    adj = 0  # for all other datasets
    
    arr_have = []
    for i in range(0,l):
        x = flist[i]
        h = int(x[(l_fname-33+adj):(l_fname-31+adj)])
        m = int(x[(l_fname-30+adj):(l_fname-28+adj)])
        s = int(x[(l_fname-27+adj):(l_fname-25+adj)])        
        t = ('%s''%02d'':''%02d'':''%02d' % (Y1,h,m,s))
        arr_have.append(t)
    #print arr_have
    print len(arr_have)    
    
    
    f_first = flist[0]
    h3 = int(f_first[(l_fname-33+adj):(l_fname-31+adj)])
    m3 = int(f_first[(l_fname-30+adj):(l_fname-28+adj)]) - 1 
    s3 = (60 +  (int(f_first[(l_fname-27+adj):(l_fname-25+adj)])) - (int(np.floor(60/cadence))*cadence))
    #t_first = ('%s''%02d'':''%02d'':''%02d' % (Y1,H1,(m_first - 1),s_first))
    
    # create array of all possible files
    arr_all = []   
    for i in range(0,num_files):  # changed from range(0,l)
        t3 = ('%s''%02d'':''%02d'':''%02d' % (Y1,h3,m3,s3))
        arr_all.append(t3)
        s3 += cadence
        if s3 >= 60:
            s3 -= 60
            m3 += 1
            if m3 >= 60:
                m3 -= 60
                h3 += 1 
    #print arr_all
    print len(arr_all)
    
    
    # compare array_all to array_have to determine array_need
    arr_need = []   
    for i in range(0,num_files):  # changed from range(0,l)
        z = arr_all[i] in arr_have    
        if z == False:
            arr_need.append(arr_all[i])
    print arr_need
    print len(arr_need)    
    
    print ""
    print "After the initial pass, still need %d files." % len(arr_need)
    
    # loop through the array of needed files, requesting them one at a time         
    for i in range(0,len(arr_need)):
        qr=client.query(vso.attrs.Time(arr_need[i],arr_need[i]), vso.attrs.Instrument('aia'), vso.attrs.Wave(wavelength * u.AA, wavelength * u.AA))
        print qr
        res=client.get(qr, path='%s/{file}.fits' % path_name).wait()
        #res=client.get(qr, path='%s/{file}' % path_name).wait()  # only for 193  (maybe for all going forward?)
        


"""
############################
############################
# generate heatmaps
############################
############################
"""

## currently enter wavelength and date as arguments in call - could possibly
## extract that info from parameter-array file name? 

# think I can take out the trimming of edges?  should be built into the array creation

# might want to flip all heatmaps / visual images to match sunpy's peek()?

# maybe instead of NaN for chi^2, use [if.data?] - so only looks at valid data

# when generating heatmaps for rebinned regions - should it scale x/y axes by that factor?

# update 1/13:
# included visual images in single function call
# added pdf-save support for easy inclusion in latex documents

# update 1/19:
# added now generates histograms of all parameters
# might want to set histogram ranges based on curve_fit parameter bounds?
# want them to be all the same ranges? -- or are we just using to pick out fit problems?

# got odd error when generating float32 heatmaps:
#UserWarning: Attempting to set identical left==right results
#in singular transformations; automatically expanding.
#left=0.425, right=0.425 'left=%s, right=%s') % (left, right))

# add in generalized - aspect ratio of region determines size of figure

# update 1/25:
# added percentiles to visual images 
# changed gaussian location ticks to seconds 

# update 1/26:
# removing plt.tight_layout gives all same size plots
# find way to possibly take aspect ratio of array and make figure size match that
# allowing for the parameters that have more digits 

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff

def heatmap(heatmaps, visual, date, wavelength, path_name):
    """
    Generates heatmaps and histograms for each of the parameters:
    Powerlaw Slope Coefficient, Powerlaw Index, Powerlaw Tail Value,
    Gaussian Amplitude, Gaussian Location, Gaussian Width, and Chi^2.
    
    heatmaps : 
        Array of parameters from which heatmaps are generated. (array)
        
    visual : 
        Array containing visual images. (array)
        
    date : 
        Date of dataset in 'YYYYMMDD' format. (String)
    
    wavelength :
        Wavelength of dataset. (Integer)
                   
    path_name : 
        The directory to which the files should be saved. (String)
      
    Example:
    ::
        ss.heatmap(heatmap = HEATMAPS, visual = VISUAL, date = '20130815',
           wavelength=211, path_name='C:/Users/Brendan/Desktop/PHYS 326') 
    """
    
    # create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
    #titles = ['Slope Coefficient', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location [sec]', 'Gaussian Width', '$\chi^2$']
    titles = [r'Power Law Slope-Coefficient -- [$A$]', r'Power Law Index -- [$n$]', r'Power Law Tail -- [$C$]', r'Gaussian Amplitude -- [$\alpha$]', r'Gaussian Location [Seconds] -- [$\beta$]', r'Gaussian Width -- [$\sigma$]', 'F-Statistic', r'Gaussian Amplitude Scaled -- [$\alpha$]', 'P-Value']
    #names = ['PL_A', 'Slopes', 'PL_C', 'Gauss_Amp', 'Gauss_Loc', 'Gauss_Wid', 'Chi2']
    names = ['PL_A', 'Slopes', 'PL_C', 'Gauss_Amp', 'Gauss_Loc', 'Gauss_Wid', 'F_Test', 'Gauss_Amp_Scaled', 'P_Value']
    #cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '$\chi^2$']
    #cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [seconds]', 'Width', '$\chi^2$']
    cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [seconds]', 'Width', 'F-Statistic', 'Amplitude Scaled', 'P-Value']
    
    #vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore  (or option to set ranges for specific wavelengths?)
    #vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore
    
    wavelength = wavelength
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    date_title = '%s-%s-%s' % (year,month,day)
    
    h_map = heatmaps
    h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]  # trim last row and column from array (originally needed since went one past)
    
    if h_map.shape[2] > h_map.shape[1]:
        aspect_ratio = float(h_map.shape[2]) / float(h_map.shape[1])
        fig_height = 10
        fig_width = 10*aspect_ratio
        
    else:
        aspect_ratio = float(h_map.shape[1]) / float(h_map.shape[2])
        #print aspect_ratio
        #fig_width = 10
        fig_width = 10+1  # works better for 20130626 (with no x/y labels)
        #fig_height = 10*aspect_ratio
        fig_height = 10*aspect_ratio  # works better for 20130626
    
    
    for i in range(0,len(titles)-1):
    #for i in range(0,1):
        
        #fig = plt.figure(figsize=(13,9))
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
        plt.title('%s' % (titles[i]), y = 1.01, fontsize=25)  # no date / wavelength
        
        if i == 6:
            NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
            h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            cmap = 'jet'
        elif i == 4:
            h_map[i] = 1./(np.exp(h_map[i]))
            h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            cmap = 'jet_r'  # reverse color-scale for Gaussian Location, because of flipped frequencies to seconds
        elif i == 8:
            df1, df2 = 3, 6
            h_map[6] = ff.sf(h_map[6], df1, df2)
            h_min = np.percentile(h_map[6],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[6],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            cmap = 'jet'                     
        else:
            h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            cmap = 'jet'
        
        im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=h_min, vmax=h_max)
        #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
        #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
        cax = divider.append_axes("right", size="3%", pad=0.07)
        if i == 0:
            cbar = plt.colorbar(im,cax=cax, format='%0.2e')
        else:
            cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=17, pad=5) 
        #plt.tight_layout()
        #plt.savefig('%s/%s_%i_heatmap_%s.jpeg' % (path_name, date, wavelength, names[i]))
        plt.savefig('%s/%s_%i_heatmap_%s.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
        
        """
        flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
    
        fig = plt.figure(figsize=(12,9))
        plt.title(r'%s: %i $\AA$  [Histogram - %s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
        plt.xlabel('%s' % cbar_labels[i], fontsize=20, labelpad=10)
        plt.ylabel('Bin Count', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        plt.xlim(h_min, h_max)
        y, x, _ = plt.hist(flat_param, bins=200, range=(h_min, h_max))
        plt.ylim(0, y.max()*1.1)
        #plt.hist(flatten_slopes, bins='auto')  # try this (actually think we want constant bins throughout wavelengths)
        #plt.savefig('%s/%s_%i_Histogram_%s.jpeg' % (path_name, date, wavelength, names[i]))
        plt.savefig('%s/%s_%i_Histogram_%s.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
        """
    
    # generate p-value heatmap
    df1, df2 = 3, 6
    p_val = ff.sf(h_map[6], df1, df2)

    p_mask = np.copy(p_val)
    
    mask_thresh = 0.005    
       
    p_mask = np.copy(p_val)
    amp_mask = np.copy(h_map[3])
    loc_mask = np.copy(h_map[4])
    wid_mask = np.copy(h_map[5])
    
    h_min_amp = np.percentile(h_map[3],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max_amp = np.percentile(h_map[3],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    
    for i in range(p_val.shape[0]):
            for j in range(p_val.shape[1]):
                if p_val[i][j] > mask_thresh:
                    p_mask[i][j] = np.NaN
                    amp_mask[i][j] = np.NaN
                    loc_mask[i][j] = np.NaN
                    wid_mask[i][j] = np.NaN
                    
    plots = [p_mask, amp_mask, loc_mask, wid_mask]
    names_f = ['P-Value Mask', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width']
    names_m = ['p_mask', 'amp', 'loc', 'wid']
    
    for k in range(4):           
    
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        if k == 0:
            plt.title('P-Value < %0.3f' % (mask_thresh), y = 1.01, fontsize=25)
        else:
            plt.title('%s: P-Value < %0.3f' % (names_f[k], mask_thresh), y = 1.01, fontsize=25)
        if k == 2:
            cmap = 'jet_r'
        else:
            cmap = 'jet'
        if k == 1:
            im = ax.imshow(np.flipud(plots[k]), cmap = cmap, vmin = h_min_amp, vmax = h_max_amp)
        else:
            im = ax.imshow(np.flipud(plots[k]), cmap = cmap)
        #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
        #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=17, pad=5) 
        #plt.tight_layout()
        #plt.savefig('%s/%s_%i_heatmap_%s.jpeg' % (path_name, date, wavelength, names[i]))
        plt.savefig('%s/%s_%i_%s_mask_%i.pdf' % (path_name, date, wavelength, names_m[k], (1./mask_thresh)), format='pdf')
        
    

  
    # generate visual images
    titles_vis = ['Average', 'Middle-File']
    names_vis = ['average', 'mid']
    
    vis = visual
    
    for i in range(2):
        
        v_min = np.percentile(vis[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        v_max = np.percentile(vis[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)     
        
        #fig = plt.figure(figsize=(12,9))
        fig = plt.figure(figsize=(fig_width,fig_height))
        
        ax = plt.gca()
        #plt.title(r'%s: %i $\AA$  [Visual: %s]' % (date_title, wavelength, titles_vis[i]), y = 1.01, fontsize=25)
        plt.title('Visual: %s' % (titles_vis[i]), y = 1.01, fontsize=25)  # no date / wavelength
        #im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
        im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
        #plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
        #plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        #cbar.set_label('Intensity', size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=17, pad=5) 
        #plt.tight_layout()
        #plt.savefig('%s/%s_%i_visual_%s.jpeg' % (path_name, date, wavelength, names_vis[i]))
        plt.savefig('%s/%s_%i_visual_%s.pdf' % (path_name, date, wavelength, names_vis[i]), format='pdf')




"""
############################
############################
# pixel-to-arcsecond + subregion check
############################
############################
"""

# doesn't account for derotation amount - build that in


import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import sunpy
from sunpy.map import Map

def arc2pix(x1, x2, y1, y2, image):
    """
    Converts arcsecond coordinates to pixel coordinates + generates subregion image.
    
    x1 : 
        Arcsecond coordinate to be converted (int/float)
        
    x2/y1/y2 : 
        Same as x1
    
    image : 
        .fits file to generate subregion image from. (String)
      
    Example:
    ::
        ss.arc2pix(x1,x2,y1,y2, image = 'F:/SDO/data/20130815/171/aia_lev1_171a_2013_08_15t05_59_59_34z_image_lev1.fits.fits') 
    """
    x1 = x1
    x2 = x2
    y1 = y1
    y2 = y2
    
    m1 = Map('%s' % image)
    x_coord = m1.data_to_pixel(x1*u.arcsec,x2*u.arcsec,1)
    x1_coord = x_coord[0].value
    x2_coord = x_coord[1].value
    y_coord = m1.data_to_pixel(y1*u.arcsec,y2*u.arcsec,1)
    y1_coord = y_coord[0].value
    y2_coord = y_coord[1].value
    m1.peek()
    m2 = Map('%s' % image).submap([x1,x2]*u.arcsec, [y1,y2]*u.arcsec)
    m2.peek()
    print "[Arcsecond Coordinate] = [Pixel Coordinate]"
    print "%i arcsec = x1 = %i pixel" % (x1,x1_coord)
    print "%i arcsec = x2 = %i pixel" % (x2,x2_coord)
    print "%i arcsec = y1 = %i pixel" % (y1,y1_coord)
    print "%i arcsec = y2 = %i pixel" % (y2,y2_coord)
    #b3 = np.flipud(b2)
    #plt.imshow(b3, cmap='sdoaia171')


def pix2arc(x1, x2, y1, y2, image):
    """
    Converts pixel coordinates to arcsecond coordinates + generates subregion image.
    
    x1 : 
        Pixel coordinate to be converted (int/float)
        
    x2/y1/y2 : 
        Same as x1
    
    image : 
        .fits file to generate subregion image from. (String)
      
    Example:
    ::
        ss.pix2arc(x1,x2,y1,y2, image = 'F:/SDO/data/20130815/171/aia_lev1_171a_2013_08_15t05_59_59_34z_image_lev1.fits.fits') 
    """
    x1 = x1
    x2 = x2
    y1 = y1
    y2 = y2
    
    m1 = Map('%s' % image)
    x_coord = m1.pixel_to_data(x1*u.pixel,x2*u.pixel,1)
    x1_coord = x_coord[0].value
    x2_coord = x_coord[1].value
    y_coord = m1.pixel_to_data(y1*u.pixel,y2*u.pixel,1)
    y1_coord = y_coord[0].value
    y2_coord = y_coord[1].value
    m1.peek()
    m2 = Map('%s' % image).submap([x1,x2]*u.pixel, [y1,y2]*u.pixel)
    m2.peek()
    print "[Pixel Coordinate] = [Arcsecond Coordinate]"
    print "%i pixel = x1 = %i arcsec" % (x1,x1_coord)
    print "%i pixel = x2 = %i arcsec" % (x2,x2_coord)
    print "%i pixel = y1 = %i arcsec" % (y1,y1_coord)
    print "%i pixel = y2 = %i arcsec" % (y2,y2_coord)    



"""
############################
############################
# datacube + derotate 
############################
############################
"""

# could get first .fits file and extract date and wavelength from
# worked for one full run-through - then also worked through full FFT + 3x3 + fit

# maybe add -- if rebin = 1, don't add to filename
# have visual / time be either both at beginning or end of filename

# add printout of region dimensions and rebinned region dimensions

from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
import numpy as np
import astropy.units as u



def datacube(directory, date, wavelength, sub_reg_coords, coords_type, bin_frac):
    """
    Generates derotated datacube of subregion.
    
    directory : 
        Path to folder containing .fits files - also where cube is saved.  (String)
        
    date :
        Date of timeseries [YYYY/MM/DD]  (String)
        
    wavelength :
        Wavelength of dataset.  (Int)
        
    sub_reg_coords : 
        Coordinates of subregion [x1,x2,y1,y2] (int/float)
    
    coords_type :
        Arcsecond or Pixel coordinates ['arc' or 'pix'] (String)
                   
    bin_frac : 
        The fraction by which the image should be rebinned by. (int)
      
    Example:
    ::
        ss.datacube(directory='F:/SDO/data/20130530/1600', date='20130530', wavelength=1600,
           sub_reg_coords=[2200,3000,2300,2600], coords_type='pix', bin_frac=2) 
    """
    
    # rebin region to desired fraction 
    def rebin(a, *args):
        shape = a.shape
        lenShape = len(shape)
        factor = np.asarray(shape)/np.asarray(args)
        evList = ['a.reshape('] + \
                 ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
                 [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
        #print ''.join(evList)
        return eval(''.join(evList))
        
    # define subregion coordinates   
    # or possibly import these from the return values of arc2pix  <-- this
    x1 = sub_reg_coords[0]  # j1
    x2 = sub_reg_coords[1]  # j2
    y1 = sub_reg_coords[2]  # i1
    y2 = sub_reg_coords[3]  # i2
    
    # create a list of all the files. This is USER-DEFINED
    flist = glob.glob('%s/aia*.fits' % directory)
    nf = len(flist)

    # Select the image that is the "middle" of our selection.
    # We do this because the solar derotation algorithm operates centered on the 
    # "middle" image  (might not be used)
    mid_file = np.int(np.floor(nf / 2))

 
    
    mc_list = []  # create an empty list
    
    count = 0  # counter to see where program is at 
    
    # Use defined coordinates, extract the submaps from each AIA image, and store
    # them in the empty list. This takes many minutes to complete.
    print " "
    print "Reading files and extracting submaps. This takes a while..."
    print " "
   
    if coords_type == 'pix':
        for filename in flist:
            mc_list.append(Map(filename).submap([x1,x2]*u.pixel, [y1,y2]*u.pixel))
            if (count%10) == 0:
                print("file %i out of %i" % (count,nf))  # just to see where program is at
            count += 1
    
    if coords_type == 'arc':
        for filename in flist:
            mc_list.append(Map(filename).submap([x1,x2]*u.arcsec, [y1,y2]*u.arcsec))
            if (count%10) == 0:
               print("file %i out of %i" % (count,nf))  # just to see where program is at
            count += 1
        
    	
    # Create a new Map Cube object to hold our de-rotated data
    new_mapcube = Map(mc_list, cube=True)
    print " "
    print "Creating derotated cube..."
    print "Please wait..."
    
    # Perform the derotation of the submaps. This take a while too.
    dr = mapcube_solar_derotate(new_mapcube)
    
    print "done derotating"
    
    
    mid_subarr = dr[mid_file].data		# extract data from middle file of derotated datacube
    rem_i = mid_subarr.shape[0] % bin_frac  # calculate remainder to get integer dimensions
    rem_j = mid_subarr.shape[1] % bin_frac  # calculate remainder to get integer dimensions
    subarr_idim = (mid_subarr.shape[0] - rem_i) / bin_frac  # get ydim
    subarr_jdim = (mid_subarr.shape[1] - rem_j) / bin_frac  # get xdim
     
    mc_list = np.zeros((1,2,3))  # free up memory
    new_mapcube = np.zeros((1,2,3))  # free up memory
    
    t = dr[0].date  # extract the date / time from the first image
    base_time = (t.hour * 60.) + (t.minute * 60.) + t.second  # convert date / time to seconds

    # initialize arrays to hold exposure time, pixel data, and time values
    I = np.empty((nf))  
    DATA = np.empty((nf, subarr_idim, subarr_jdim))
    #TIME = np.empty((nf), dtype='int')
    TIME = np.empty((nf))  # might as well just have all as float
    
    # loop through datacube and extract pixel data and time values
    for p in range(0,nf):
        Ex = dr[p].exposure_time
        I[p] = Ex.value
        L = dr[p].data
        L_trim = L[0:(mid_subarr.shape[0] - rem_i), 0:(mid_subarr.shape[1] - rem_j)]
        small_L = rebin(L_trim, L_trim.shape[0]/bin_frac, L_trim.shape[1]/bin_frac)
        DATA[p][:][:] = small_L/Ex  # normalize by exposure time
        T = dr[p].date
        curr_time=(T.hour * 3600.)+(T.minute * 60.)+T.second	
        TIME[p] = curr_time - base_time  # calculate running time of image
    
    # save the pixel-value and time-array datacubes as numpy files
    np.save('%s/%s_%i_%i_%ii_%i_%ij_data_rebin%i.npy' % (directory, date, wavelength, y1, y2, x1, x2, bin_frac), DATA)
    np.save('%s/%s_%i_%i_%ii_%i_%ij_time.npy' % (directory, date, wavelength, y1, y2, x1, x2), TIME)
    
    # calculate the average-intensity image of the timeseries 
    AVG = np.average(DATA,axis=0)
    
    # determine the middle file of the timeseries
    mid_num = (DATA.shape[0]/2)
    mid = DATA[mid_num]
    
    print "Middle file is number %i" % mid_num
    
    print "(100,100): %i ~ %i" % (AVG[100][100], mid[100][100])  # check values are reasonably close
    
    # check the two image sizes agree
    print " Average Image Dimensions = %i, %i" % (AVG.shape[0], AVG.shape[1])
    print " Middle Image Dimensions = %i, %i" % (mid.shape[0], mid.shape[1])
    
    # store average and middle images in array
    visual = np.zeros((2,AVG.shape[0],AVG.shape[1]))
    visual[0] = AVG
    visual[1] = mid
    
    print visual.shape  # check array size agrees with expected
    
    # generate images of each visual region, to see if as expected
    fig = plt.figure(figsize=(20,20))
    plt.imshow(visual[0])
    fig = plt.figure(figsize=(20,20))
    plt.imshow(visual[1])
    
    # save visual-image array
    np.save('%s/visual_%s_%i_%i_%ii_%i_%ij.npy' % (directory, date, wavelength, y1, y2, x1, x2), visual)



"""
############################
############################
# FFT segment averaging + 3x3 Pixel Box Averaging
############################
############################
"""

# update 1/16: including 3x3 pixel-box averaging within this function

# update 1/18: using accelerate, FFT speedup by 2x!


import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate  # switch on if computer has installed


def fft_avg(datacube, timeseries, num_seg):
    """
    Calculates segment-averaged FFT for region, and 3x3 pixel-box average.
    
    datacube : 
        3D datacube of timeseries of images.  (array)
        
    timeseries :
        Array of timeseries.  (array)
        
    num_seg :
        Number of segments to divide timeseries into.  (int)
        
    Returns : 
        3D array of FFT-segment and 3x3-pixel-box averaged region
      
    Example:
    ::
        ss.fft_avg(datacube = DATA, timeseries = TIME, num_seg = 3)
    """
    
    # currently not using
    #with h5py.File('C:/Users/Brendan/Desktop/SDO/SDO_20120923_211A_(528)_(132)x_(100)_100y.hdf5', 'r', driver='core') as f:
        #DATA = f['Image Data']
        #TIME = f['Time Data']
    
    from scipy import fftpack
    
    DATA = datacube
    
    TIME = timeseries
    
    print DATA.shape 
    
    print "Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0])
    
    t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters
    
    n_segments = num_seg  # break data into 12 segments of equal length
    n = len(t_interp)
    rem = n % n_segments
    freq_size = (n - rem) / n_segments
    
    ## determine frequency values that FFT will evaluate
    time_step = 12  # add as argument in function call, or leave in as constant?
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    reslt = (DATA.shape[0] == TIME.shape[0])
    print "DATA and TIME array sizes match: %s" % reslt
    
    pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
    spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))
    
    print "length time-interp array = %i" % n
    print "size for FFT to consider = %i" % freq_size
    print "length of sample freq array = %i" % len(sample_freq)
    print "length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs)
    
    
    start = timer()
    T1 = 0
    
    for ii in range(0,spectra_seg.shape[0]):
    #for ii in range(142,160):
    
        for jj in range(0,spectra_seg.shape[1]):
        #for jj in range(0,574):        
        
            x1_box = 0+ii
            #x2_box = 2+ii  # if want to use median of more than 1x1 pixel box
            y1_box = 0+jj
            #y2_box = 2+jj  # if want to use median of more than 1x1 pixel box
            
            for k in range(0,DATA.shape[0]):
                        im=DATA[k]	  # get image
                        #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])  # finds pixel-box median
                        pixmed[k]= im[x1_box,y1_box]	# median  <-- use this
                        
                        
            # The derotation introduces some bad data towards the end of the sequence. This trims that off
            bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
            last_good_pos = bad - 1			# retain only data before the <=zero
            
            # Get time and pixel values
            v=pixmed[0:last_good_pos]		
            t=TIME[0:last_good_pos]
            
            #plt.plot(t,v)
        
            v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
            
            data = v_interp
            
            avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
    
            data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
            split = np.split(data, n_segments)  # create split array for each segment
    
    
            for i in range(0,n_segments):               
                
                 ## perform Fast Fourier Transform on each segment       
                 sig = split[i]
                 sig_fft = fftpack.fft(sig)
                 #sig_fft = fftpack.rfft(sig)  # real-FFT
                 #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
                 #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
                 #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
                 powers = np.abs(sig_fft)[pidxs]
                 norm = len(sig)  # to normalize the power
                 powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
                 avg_array += powers
            
            avg_array /= n_segments  # take the average of the segments
            
            avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
                       
            spectra_seg[ii][jj] = avg_array  # construct 3D array with averaged FFTs from each pixel
        
        
        # estimate time remaining and print to screen
        T = timer()
        T2 = T - T1
        if ii == 0:
            T_init = T - start
            T_est = T_init*(spectra_seg.shape[0])  
            print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est)
        else:
            T_est2 = T2*(spectra_seg.shape[0]-ii)
            print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est2)
        T1 = T
        
        
    # print estimated and total program time to screen 
    print "Beginning Estimated time = %i sec" % T_est
    T_act = timer() - start
    print "Actual total time = %i sec" % T_act 
    
    
    # initialize arrays to hold temporary results for calculating geometric average
    temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
    p_geometric = np.zeros((spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
    spectra_array = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
        
    
    ### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
    ## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

    for l in range(1,spectra_seg.shape[0]-1):
    #for l in range(1,25):
        #print l
        for m in range(1,spectra_seg.shape[1]-1):
        #for m in range(1,25):
            
            temp[0] = np.log10(spectra_seg[l-1][m-1])
            temp[1] = np.log10(spectra_seg[l-1][m])
            temp[2] = np.log10(spectra_seg[l-1][m+1])
            temp[3] = np.log10(spectra_seg[l][m-1])
            temp[4] = np.log10(spectra_seg[l][m])
            temp[5] = np.log10(spectra_seg[l][m+1])
            temp[6] = np.log10(spectra_seg[l+1][m-1])
            temp[7] = np.log10(spectra_seg[l+1][m])
            temp[8] = np.log10(spectra_seg[l+1][m+1])
    
            temp9 = np.sum(temp, axis=0)
            p_geometric = temp9 / 9.
            spectra_array[l-1][m-1] = np.power(10,p_geometric)
    
    return spectra_array
    


"""
############################
############################
# Curve Fitting
############################
############################
"""    

# update 1/14:
# changing slope bounds to [0.3, 4.0]
# moved the 3x3 pixel-box averaging, only fitting now

# when segmenting for parallelization - numpy.load(mmap-reads in only slice?)

# should revise the variable naming ('gp' --> 'M2', all the '2' variables?)

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import h5py
from scipy import fftpack  # doesnt work in module when called here???
import matplotlib.pylab as plt
from astropy.convolution import convolve, Box1DKernel
from matplotlib import cm
from numpy.random import randn
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.ticker import LogFormatterMathtext
from timeit import default_timer as timer
from scipy.stats import f


def spec_fit(spectra_array):
    """
    Calculates pixel-box averages and fits spectra.
    
    spectra_array : 
        3D array of FFT-segment and 3x3-pixel-box averaged region  (array)
        
    Returns : 
        Parameter and M2-fit arrays
      
    Example:
    ::
        ss.spec_fit(spectra_array = SPECTRA)
    """

    from scipy import fftpack    
    
    # define Power-Law-fitting function (Model M1)
    def PowerLaw(f, A, n, C):
        return A*f**-n + C
        
    # define Gaussian-fitting function
    def Gauss(f, P, fp, fw):
        return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)
    
    # define combined-fitting function (Model M2)
    def GaussPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
        return A2*f2**-n2 + C2 + P2*np.exp(-0.5*(((np.log(f2))-fp2)/fw2)**2)
        
    
    ## load in array of segment-averaged pixel FFTs
    SPECTRA = spectra_array
    #spectra_array = spectra_array.astype(np.float32)  # possibly use this?
    
    print "The region size is %ii x %ij" % (SPECTRA.shape[0], SPECTRA.shape[1])
    print "%i frequencies were evaluated in the FFT" % SPECTRA.shape[2] 
    
    num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used
        
    # determine frequency values that FFT will evaluate
    freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
    time_step = 12  # add as argument, or leave in as constant?
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    #freqs = freqs.astype(np.float32)  # possibly use this?
    
    
    # initialize arrays to hold parameter values, also each pixel's combined model fit - for tool
    #diffM1M2 = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1]))  # dont really use - get rid of?
    params = np.zeros((8, SPECTRA.shape[0], SPECTRA.shape[1]))
    #M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], (len(freqs)+1)/2))  # would save storage / memory space
    M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], SPECTRA.shape[2]))
    
    #Uncertainties = np.zeros((6, SPECTRA.shape[0], SPECTRA.shape[1]))
    
    start = timer()
    T1 = 0
    
    
    ### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
    ## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)
    
    for l in range(0,SPECTRA.shape[0]):
    #for l in range(145,146):
        
        for m in range(0,SPECTRA.shape[1]):
        #for m in range(185,195):
            
                                            
            f = freqs  # frequencies
            s = spectra_array[l][m]  # fourier power
            
            #ds = (1./f**2.2)/1000
            ds = s*0.1  # set the error / variance estimate to a constant percentage of the spectra power-values
            
            # create points to fit model with final parameters 
            #f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],(len(freqs)+1)/2)  # would save storage / memory space?
            #f_fit = freqs       
            
                                                   
            ### fit data to models using SciPy's Levenberg-Marquart method
            
            try:
                # initial guesses for fitting parameters
                M1_low = [-0.002, 0.3, -0.01]
                M1_high = [0.002, 4., 0.01]
                nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')  # replaced #'s with arrays
               
            
            except RuntimeError:
                print("Error M1 - curve_fit failed - %i, %i" % (l,m))
            
            except ValueError:
                print("Error M1 - inf/NaN - %i, %i" % (l,m))
    
          
            A, n, C = nlfit_l  # unpack fitting parameters
            
            # unpack uncertainties in fitting parameters from diagonal of covariance matrix
            dA, dn, dC = [np.sqrt(nlpcov_l[j,j]) for j in range(nlfit_l.size)]
            
            
            ## possibly use previous pixel's parameters as initial guesses for current pixel (issues creating wierd banding in images)
                            
            """
            M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
            #M2_low = [-0.1, 0.1, -0.1, 0.00001, -6.5, 0.05]  # test on 193 - coronal hole
            M2_high = [0.002, 4., 0.01, 0.2, -4.6, 0.8]
            if m > 0:
                P0 = [params[0][l][m-1], params[1][l][m-1], params[2][l][m-1], params[3][l][m-1], params[4][l][m-1], params[5][l][m-1]]
                #P0 = [0.0, 1.02, 0.001, 0.001, -4.68, 0.79]
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
                    
            elif m == 0 and l > 0:
                P0 = [params[0][l-1][m], params[1][l-1][m], params[2][l-1][m], params[3][l-1][m], params[4][l-1][m], params[5][l-1][m]]
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
                    
            else:        
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0=[nlfit_l[0], nlfit_l[1], nlfit_l[2], nlfit_g[0], nlfit_g[1], nlfit_g[2]], sigma=ds)
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
            """     
            
            ## fit data to combined power law plus gaussian component model
            #"""        
            try:                                 
                M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
                M2_high = [0.002, 4., 0.01, 0.2, -4.6, 0.8]
                #M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]  # see what happens if force middle of range above where slopes are
                
                # change method to 'dogbox' and increase max number of function evaluations to 3000
                nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000) # replaced #'s with arrays
                
            except RuntimeError:
                print("Error M2 - curve_fit failed - %i, %i" % (l,m))
            
            except ValueError:
                print("Error M2 - inf/NaN - %i, %i" % (l,m))
            #"""
            
            A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
            
            # unpack uncertainties in fitting parameters from diagonal of covariance matrix
            dA2, dn2, dC2, dP2, dfp2, dfw2 = [np.sqrt(nlpcov_gp[j,j]) for j in range(nlfit_gp.size)]
            
            m2_param = A2, n2, C2, P2, fp2, fw2  # could have used this for params array : = params[0:6,l-1,m-1]
            uncertainties = dA2, dn2, dC2, dP2, dfp2, dfw2  # do we want to keep a global array of uncertainties?  
            
            
            #uncertainties_arr = [dA2, dn2, dC2, dP2, dfp2, dfw2]  # not sure if want to keep these
            #Uncertainties[:, l, m] = uncertainties_arr
            
            
            # create model functions from fitted parameters
            #m1_fit = PowerLaw(f_fit, A, n, C)
            m1_fit = PowerLaw(f, A, n, C)
            amp_scale = PowerLaw(np.exp(fp2), A, n, C)  # to extract the gaussian-amplitude scaling factor
            #m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
            m2_fit = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
            #s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)  # could get rid of this if not making smaller m2_fit
            #m2P_fit = PowerLaw(f_fit, A2, n2, C2)  # only need if plotting
            #m2G_fit = Gauss(f_fit, P2, fp2, fw2)  # only need if plotting
            
            #diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
            #diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
                                   
            
            #residsM2 = (s - s_fit_gp_full)
            residsM2 = (s - m2_fit)
            chisqrM2 = ((residsM2/ds)**2).sum()
            redchisqrM2 = ((residsM2/ds)**2).sum()/float(f.size-6)
            
            residsM1 = (s - m1_fit)
            chisqrM1 =  ((residsM1/ds)**2).sum()
            redchisqrM1 = ((residsM1/ds)**2).sum()/float(f.size-3)       
            
            f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))
            
            # populate array with parameters
            params[0][l][m] = A2
            params[1][l][m] = n2
            params[2][l][m] = C2
            params[3][l][m] = P2
            params[4][l][m] = fp2
            params[5][l][m] = fw2
            #params[6][l][m] = redchisqrM2
            params[6][l][m] = f_test
            params[7][l][m] = P2 / amp_scale
            
            
            # populate array holding model fits
            M2_fit[l][m] = m2_fit
            
            
            
            # Plot models + display combined-model parameters + uncertainties
            """
            fig = plt.figure(figsize=(20,15))
            plt.title('SDO AIA 304.0 Angstrom 20120923 - 3 Segments, 3x3 Pixel-Box Averaging 598 interp', y = 1.01, fontsize=25)
            plt.ylim((10**-8,10**1))
            plt.xlim((10**-5,10**-1))
            plt.loglog(f,s,'k')
            plt.loglog(f_fit, m1_fit, label='Power Law - M1')
            plt.loglog(f_fit, m2P_fit, 'g', label='Power Law - M2')
            plt.loglog(f_fit, m2G_fit, 'g--', label='Gaussian - M2')
            plt.loglog(f_fit, m2_fit, 'r', label='Combined - M2')
            plt.xlabel('Frequency (Hz)', fontsize=20, labelpad=10)
            plt.ylabel('Power', fontsize=20, labelpad=10)
            plt.vlines((1.0/300.),10**-8,10**1, linestyles='dashed', label='5 minutes')
            plt.vlines((1.0/180.),10**-8,10**1, linestyles='dotted', label='3 minutes')
            plt.vlines((1.0/24.),10**-8,10**1, linestyles='solid', label='24 seconds')
            plt.text(0.01, 10**0., 'A = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[0], uncertainties[0]), fontsize=15)
            plt.text(0.01, 10**-0.3, 'n = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[1], uncertainties[1]), fontsize=15)
            plt.text(0.01, 10**-0.6, 'C = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[2], uncertainties[2]), fontsize=15)
            plt.text(0.01, 10**-0.9, 'P = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[3], uncertainties[3]), fontsize=15)
            plt.text(0.01, 10**-1.2, 'fp = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[4], uncertainties[4]), fontsize=15)
            plt.text(0.01, 10**-1.5, 'fw = {0:0.3f}$\pm${1:0.3f}'.format(m2_param[5], uncertainties[5]), fontsize=15)
            plt.text(0.01, 10**-1.8, 'f_test = {0:0.3f}'.format(f_test), fontsize=15)
            plt.legend(loc='upper left', prop={'size':15})
            plt.show()
            #plt.savefig('C:/Users/Brendan/Desktop/PHYS 326/dogbox_test1/20130530_193A_3x3_6seg_%ii_%ij.jpeg' % (l,m))
            #plt.savefig('C:/Users/Brendan/Desktop/SDO/20120923_%ii_%ij_598_interp.jpeg' % (l,m))
            #plt.close()
            """
        
        # estimate time remaining and print to screen  (looks to be much better - not sure why had above?)
        T = timer()
        T2 = T - T1
        if l == 0:
            T_init = T - start
            T_est = T_init*(SPECTRA.shape[0])  
            print "Currently on row %i of %i, estimated time remaining: %i seconds" % (l, SPECTRA.shape[0], T_est)
        else:
            T_est2 = T2*((SPECTRA.shape[0])-l)
            print "Currently on row %i of %i, estimated time remaining: %i seconds" % (l, SPECTRA.shape[0], T_est2)
        T1 = T
    
    # print estimated and total program time to screen        
    print "Beginning Estimated time = %i sec" % T_est
    T_act = timer() - start
    print "Actual total time = %i sec" % T_act           
    
    return params, M2_fit
    
    
    
"""
############################
############################
# datacube + derotate (int) 
############################
############################
"""

# could get first .fits file and extract date and wavelength from
# worked for one full run-through - then also worked through full FFT + 3x3 + fit

# maybe add -- if rebin = 1, don't add to filename
# have visual / time be either both at beginning or end of filename

# add printout of region dimensions and rebinned region dimensions

from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
import numpy as np
import astropy.units as u



def datacube_int(directory, date, wavelength, sub_reg_coords, coords_type, bin_frac):
    """
    Generates derotated datacube of subregion.
    
    directory : 
        Path to folder containing .fits files - also where cube is saved.  (String)
        
    date :
        Date of timeseries [YYYY/MM/DD]  (String)
        
    wavelength :
        Wavelength of dataset.  (Int)
        
    sub_reg_coords : 
        Coordinates of subregion [x1,x2,y1,y2] (int/float)
    
    coords_type :
        Arcsecond or Pixel coordinates ['arc' or 'pix'] (String)
                   
    bin_frac : 
        The fraction by which the image should be rebinned by. (int)
      
    Example:
    ::
        ss.datacube(directory='F:/SDO/data/20130530/1600', date='20130530', wavelength=1600,
           sub_reg_coords=[2200,3000,2300,2600], coords_type='pix', bin_frac=2) 
    """
    
    # rebin region to desired fraction 
    def rebin(a, *args):
        shape = a.shape
        lenShape = len(shape)
        factor = np.asarray(shape)/np.asarray(args)
        evList = ['a.reshape('] + \
                 ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
                 [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
        #print ''.join(evList)
        return eval(''.join(evList))
        
    # define subregion coordinates   
    # or possibly import these from the return values of arc2pix  <-- this
    x1 = sub_reg_coords[0]  # j1
    x2 = sub_reg_coords[1]  # j2
    y1 = sub_reg_coords[2]  # i1
    y2 = sub_reg_coords[3]  # i2
    
    # create a list of all the files. This is USER-DEFINED
    flist = glob.glob('%s/aia*.fits' % directory)
    nf = len(flist)

    # Select the image that is the "middle" of our selection.
    # We do this because the solar derotation algorithm operates centered on the 
    # "middle" image  (might not be used)
    mid_file = np.int(np.floor(nf / 2))

 
    
    mc_list = []  # create an empty list
    
    count = 0  # counter to see where program is at 
    
    # Use defined coordinates, extract the submaps from each AIA image, and store
    # them in the empty list. This takes many minutes to complete.
    print " "
    print "Reading files and extracting submaps. This takes a while..."
    print " "
   
    if coords_type == 'pix':
        for filename in flist:
            mc_list.append(Map(filename).submap([x1,x2]*u.pixel, [y1,y2]*u.pixel))
            if (count%10) == 0:
                print("file %i out of %i" % (count,nf))  # just to see where program is at
            count += 1
    
    if coords_type == 'arc':
        for filename in flist:
            mc_list.append(Map(filename).submap([x1,x2]*u.arcsec, [y1,y2]*u.arcsec))
            if (count%10) == 0:
               print("file %i out of %i" % (count,nf))  # just to see where program is at
            count += 1
        
    	
    # Create a new Map Cube object to hold our de-rotated data
    new_mapcube = Map(mc_list, cube=True)
    print " "
    print "Creating derotated cube..."
    print "Please wait..."
    
    # Perform the derotation of the submaps. This take a while too.
    dr = mapcube_solar_derotate(new_mapcube)
    
    print "done derotating"
    
    
    mid_subarr = dr[mid_file].data		# extract data from middle file of derotated datacube
    rem_i = mid_subarr.shape[0] % bin_frac  # calculate remainder to get integer dimensions
    rem_j = mid_subarr.shape[1] % bin_frac  # calculate remainder to get integer dimensions
    subarr_idim = (mid_subarr.shape[0] - rem_i) / bin_frac  # get ydim
    subarr_jdim = (mid_subarr.shape[1] - rem_j) / bin_frac  # get xdim
     
    mc_list = np.zeros((1,2,3))  # free up memory
    new_mapcube = np.zeros((1,2,3))  # free up memory
    
    t = dr[0].date  # extract the date / time from the first image
    base_time = (t.hour * 60.) + (t.minute * 60.) + t.second  # convert date / time to seconds

    # initialize arrays to hold exposure time, pixel data, and time values
    I = np.empty((nf))  # exposure time
    DATA = np.empty((nf, subarr_idim, subarr_jdim), dtype=np.int16)  # save as int16, since that is what original is
    #TIME = np.empty((nf), dtype='int')
    TIME = np.empty((nf))  # might as well just have all as float
    
    # loop through datacube and extract pixel data and time values
    for p in range(0,nf):
        Ex = dr[p].exposure_time
        I[p] = Ex.value
        L = dr[p].data
        L_trim = L[0:(mid_subarr.shape[0] - rem_i), 0:(mid_subarr.shape[1] - rem_j)]
        small_L = rebin(L_trim, L_trim.shape[0]/bin_frac, L_trim.shape[1]/bin_frac)
        #DATA[p][:][:] = small_L/Ex  # normalize by exposure time
        DATA[p][:][:] = small_L  # normalize by exposure time
        T = dr[p].date
        curr_time=(T.hour * 3600.)+(T.minute * 60.)+T.second	
        TIME[p] = curr_time - base_time  # calculate running time of image
    
    # save the pixel-value, time-array, and exposure-time datacubes as numpy files
    np.save('%s/%s_%i_%i_%ii_%i_%ij_data_rebin%i.npy' % (directory, date, wavelength, y1, y2, x1, x2, bin_frac), DATA)
    np.save('%s/%s_%i_%i_%ii_%i_%ij_time.npy' % (directory, date, wavelength, y1, y2, x1, x2), TIME)
    np.save('%s/%s_%i_%i_%ii_%i_%ij_exposure.npy' % (directory, date, wavelength, y1, y2, x1, x2), I)
    
    # calculate the average-intensity image of the timeseries 
    AVG = np.average(DATA,axis=0)
    
    # determine the middle file of the timeseries
    mid_num = (DATA.shape[0]/2)
    mid = DATA[mid_num]
    
    print "Middle file is number %i" % mid_num
    
    print "(100,100): %i ~ %i" % (AVG[100][100], mid[100][100])  # check values are reasonably close
    
    # check the two image sizes agree
    print " Average Image Dimensions = %i, %i" % (AVG.shape[0], AVG.shape[1])
    print " Middle Image Dimensions = %i, %i" % (mid.shape[0], mid.shape[1])
    
    # store average and middle images in array
    visual = np.zeros((2,AVG.shape[0],AVG.shape[1]))
    visual[0] = AVG
    visual[1] = mid
    
    print visual.shape  # check array size agrees with expected
    
    # generate images of each visual region, to see if as expected
    fig = plt.figure(figsize=(20,20))
    plt.imshow(visual[0])
    fig = plt.figure(figsize=(20,20))
    plt.imshow(visual[1])
    
    # save visual-image array
    np.save('%s/%s_%i_%i_%ii_%i_%ij_visual.npy' % (directory, date, wavelength, y1, y2, x1, x2), visual)
    
    
    

"""
############################
############################
# FFT segment averaging + 3x3 Pixel Box Averaging (int)
############################
############################
"""

# added 1/29:
# if choose to save datacube as type-int, plus exposure array
# this normalizes by exposure time when extracting pixel-intensity values
# should reduce datacube array size by 4x, and not cost much time here

# including exposure normalization in loop cause time to go nuts??


import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate  # switch on if computer has installed


def fft_avg_int(datacube, timeseries, exposure_array, num_seg):
    """
    Calculates segment-averaged FFT for region, and 3x3 pixel-box average.
    
    datacube : 
        3D datacube of timeseries of images.  (array)
        
    timeseries :
        Array of timeseries.  (array)
                
    exposure_array :
        Array of exposure durations.  (array)
        
    num_seg :
        Number of segments to divide timeseries into.  (int)
        
    Returns : 
        3D array of FFT-segment and 3x3-pixel-box averaged region
      
    Example:
    ::
        ss.fft_avg(datacube = DATA, timeseries = TIME, exposure_array = Ex, num_seg = 3)
    """
    
    # currently not using
    #with h5py.File('C:/Users/Brendan/Desktop/SDO/SDO_20120923_211A_(528)_(132)x_(100)_100y.hdf5', 'r', driver='core') as f:
        #DATA = f['Image Data']
        #TIME = f['Time Data']
    
    from scipy import fftpack
    
    DATA = datacube
    
    TIME = timeseries
    
    Ex = exposure_array
    
    print DATA.shape 
    
    print "Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0])
    
    t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/12)+1)  #  <-- use this (might be correct method) - not sure if matters
    
    n_segments = num_seg  # break data into 12 segments of equal length
    n = len(t_interp)
    rem = n % n_segments
    freq_size = (n - rem) / n_segments
    
    ## determine frequency values that FFT will evaluate
    time_step = 12  # add as argument in function call, or leave in as constant?
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    reslt = (DATA.shape[0] == TIME.shape[0])
    print "DATA and TIME array sizes match: %s" % reslt
    
    pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
    spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))
    
    print "length time-interp array = %i" % n
    print "size for FFT to consider = %i" % freq_size
    print "length of sample freq array = %i" % len(sample_freq)
    print "length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs)
    
    
    start = timer()
    T1 = 0
    
    for ii in range(0,spectra_seg.shape[0]):
    #for ii in range(0,5):
    
        for jj in range(0,spectra_seg.shape[1]):
        #for jj in range(0,5):        
        
            x1_box = 0+ii
            #x2_box = 2+ii  # if want to use median of more than 1x1 pixel box
            y1_box = 0+jj
            #y2_box = 2+jj  # if want to use median of more than 1x1 pixel box
            
            for k in range(0,DATA.shape[0]):
              #im=DATA[k]/(Ex[k])	  # get image + normalize by exposure time  (time went nuts?)
              im=DATA[k]
              #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])  # finds pixel-box median
              pixmed[k]= im[x1_box,y1_box]	# median  <-- use this

            pixmed = pixmed/Ex  # normalize by exposure time    
            
            # The derotation introduces some bad data towards the end of the sequence. This trims that off
            bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
            last_good_pos = bad - 1			# retain only data before the <=zero
            
            # Get time and pixel values
            v=pixmed[0:last_good_pos]		
            t=TIME[0:last_good_pos]
            #v=pixmed  # use for 335/131/094 -- can't get rid of negative values for those
            #t=TIME
            
            #plt.plot(t,v)
        
            v_interp = np.interp(t_interp,t,v)  # interpolate pixel-intensity values onto specified time grid
            
            data = v_interp
            
            avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
    
            data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
            split = np.split(data, n_segments)  # create split array for each segment
    
    
            for i in range(0,n_segments):               
                
              ## perform Fast Fourier Transform on each segment       
              sig = split[i]
              sig_fft = fftpack.fft(sig)
              #sig_fft = fftpack.rfft(sig)  # real-FFT
              #sig_fft = np.fft.rfft(sig)  # numpy significantly slower than scipy                 
              #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
              #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
              powers = np.abs(sig_fft)[pidxs]
              norm = len(sig)  # to normalize the power
              powers = ((powers/norm)**2)*(1./(sig.std()**2))*2
              avg_array += powers
            
            avg_array /= n_segments  # take the average of the segments
            
            avg_array = np.transpose(avg_array)  # take transpose of array to fit more cleanly in 3D array
                       
            spectra_seg[ii][jj] = avg_array  # construct 3D array with averaged FFTs from each pixel
        
        
        # estimate time remaining and print to screen
        T = timer()
        T2 = T - T1
        if ii == 0:
            T_init = T - start
            T_est = T_init*(spectra_seg.shape[0])  
            T_min, T_sec = divmod(T_est, 60)
            T_hr, T_min = divmod(T_min, 60)
            #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est)
            print "Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spectra_seg.shape[0], T_hr, T_min, T_sec)
        else:
            T_est2 = T2*(spectra_seg.shape[0]-ii)
            T_min2, T_sec2 = divmod(T_est2, 60)
            T_hr2, T_min2 = divmod(T_min2, 60)
            #print "Currently on row %i of %i, estimated time remaining: %i seconds" % (ii, spectra_seg.shape[0], T_est2)
            print "Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spectra_seg.shape[0], T_hr2, T_min2, T_sec2)
        T1 = T
        
        
    # print estimated and total program time to screen 
    print "Beginning Estimated time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec)
    T_act = timer() - start
    T_min3, T_sec3 = divmod(T_act, 60)
    T_hr3, T_min3 = divmod(T_min3, 60)
    print "Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3) 
    
    
    # initialize arrays to hold temporary results for calculating geometric average
    temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
    p_geometric = np.zeros((spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
    spectra_array = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))  # would pre-allocating help? (seems to)
        
    
    ### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
    ## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

    for l in range(1,spectra_seg.shape[0]-1):
    #for l in range(1,25):
        #print l
        for m in range(1,spectra_seg.shape[1]-1):
        #for m in range(1,25):
            
            temp[0] = np.log10(spectra_seg[l-1][m-1])
            temp[1] = np.log10(spectra_seg[l-1][m])
            temp[2] = np.log10(spectra_seg[l-1][m+1])
            temp[3] = np.log10(spectra_seg[l][m-1])
            temp[4] = np.log10(spectra_seg[l][m])
            temp[5] = np.log10(spectra_seg[l][m+1])
            temp[6] = np.log10(spectra_seg[l+1][m-1])
            temp[7] = np.log10(spectra_seg[l+1][m])
            temp[8] = np.log10(spectra_seg[l+1][m+1])
    
            temp9 = np.sum(temp, axis=0)
            p_geometric = temp9 / 9.
            spectra_array[l-1][m-1] = np.power(10,p_geometric)
    
    return spectra_array