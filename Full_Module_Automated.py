# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 22:30:43 2016

@author: Brendan
"""

### Full module containing all parts from downloading of data to generating of heatmaps

# still need to : 
# add all modules once generalized 
# look for optimization 
# get rid redundant module imports
# document


# so far have tested all and work

# maybe make separate function with whole program as one from get data to heatmaps?

# maybe put module imports inside function definitions (so not executing all for each)

# update 1/14: deleted visual function since had put in heatmaps

# maybe write script to compare results easier

"""
############################
############################
# download data 
############################
############################
"""

### module for downloading VSO data - specify wavelength and input / output dates

## Status : working
## need to document and comment?

# where do I specify path to save (within function or before/during/after function call)
# put in error handling if start date is after end date? also if data is entered incorrectly - not yyyy/mm/dd format
#
# generalize for 1600 - string bounds wont work - add 1

# maybe make first past be much quicker - like 10 min at a time

# for 20130815 193 - files are just .fits not .fits.fits - screws up call obviously

from sunpy.net import vso
import astropy.units as u
import glob
import numpy as np

def get_data(wavelength, time_begin, time_end, path_name):
    """
    Downloads .fits image files from database. 
    
    Parameters
    ----------
    
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
        fm.get_data(wavelength=1600, time_begin='2016/09/23 00:00:00', time_end='2016/09/23 00:05:00', 
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
    
    #l_time = 
    
    #fname_trim = len()
   
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
    
    
    arr_all = []
    
    # create array of all possible files
    #for i in range(0,num_files):  # should change to this?
    for i in range(0,l):
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
    
    arr_need = []
    
    # compare array of 'have' to array of 'all' to determine array of 'need'
    #for i in range(0,num_files):
    for i in range(0,l):
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
# determined that needed around 120 files to fill
# stopped / got stuck after getting about 20 of them - not sure why
        
# tried again later and worked better - odd, not recognizing need 12 more from 5-6 hr?
        
# update 1/10 : changed iteration length from (0,l) to (0,num_files) - should have been initially
        
# update 1/11 : added some comments
        
# prof weigel mentioned that could instead use push/pop to move first element of array to bottom
# so that program doesn't get stuck on one file
# he also mentioned that it would be better to use the ISO module for date/time
        
""" 
when redownloading - takes first and subtracts one - but if first had already been 
downloaded - screws things up
"""

        
from sunpy.net import vso
import astropy.units as u
import glob
import numpy as np

def get_data_fill(wavelength, cadence, time_begin, time_end, path_name):
    """
    Downloads .fits image files from database. 
    
    Parameters
    ----------
    
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
        fm.get_data_fill(wavelength=171, cadence=12, time_begin='2013/08/15 00:00:00',
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
    
    #l_time = 
    
    #fname_trim = len()
   
    # find first file after 00:00:00 - set time base to that
    
    # loop through flist couple times until get all.
    
    arr_have = []
    
    #adj = 5  # adjust 5 characters for 20130815 193 dataset (doesn't have extra '.fits')
    adj = 0  # for all other datasets
    
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
    
    
    arr_all = []
    
    for i in range(0,num_files):
    #for i in range(0,l):
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
    
    arr_need = []
    
    for i in range(0,num_files):
    #for i in range(0,l):
        z = arr_all[i] in arr_have    
        if z == False:
            arr_need.append(arr_all[i])
    print arr_need
    print len(arr_need)    
    
    print ""
    print "After the initial pass, still need %d files." % len(arr_need)
            
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
## have so extracted that info from parameter-array file name?
## what should bin bounds of histogram be?  

# update 12/27 - generalized size of colorbar to match size / aspect ratio of axes

# update 1/8 - generalized y-bounds on histogram

# think I can take out the trimming of edges?  should be built into the array creation

# might want to flip all heatmaps / visual images to match peek()?

# maybe instead of NaN for chi^2, use [if.data?] - only look at valid data

# maybe no interpolation as well for imshow()

# when generating heatmaps for rebinned regions - should scale x/y axes by that factor?

# update 1/13:
# included visual images in single function call
# added pdf-save support for easy inclusion in latex documents

# update 1/19:
# added now generates histograms of all parameters
# might want to set histogram range based on parameter bounds in curve_fit?
# want them to be all the same ranges? -- or just using to pick out problems?

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

def heatmap(heatmaps, visual, date, wavelength, path_name):
    """
    Generates heatmaps for each of the parameters:
    Powerlaw Slope Coefficient, Powerlaw Index, Powerlaw Tail Value,
    Gaussian Amplitude, Gaussian Location, Gaussian Width, and Chi^2.
    Also generates a histogram of the slopes in the region.
    
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
        fm.heatmap(heatmap = HEATMAPS, visual = VISUAL, date = '20130815', wavelength=211, path_name='C:/Users/Brendan/Desktop/PHYS 326') 
    """

    titles = ['Power Law Slope Coefficient', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width', '($\chi^2$)']
    names = ['PL_A', 'Slopes', 'PL_C', 'Gauss_Amp', 'Gauss_Loc', 'Gauss_Wid', 'Chi2']
    cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '($\chi^2$)']
    #vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore
    #vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore
    wavelength = wavelength
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    date_title = '%s-%s-%s' % (year,month,day)
    
    h_map = heatmaps
    h_map = h_map[:,0:h_map.shape[1]-1,0:h_map.shape[2]-1]
    
    
    for i in range(0,len(titles)):
        
        fig = plt.figure(figsize=(15,9))
        ax = plt.gca()
        plt.title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[i]), y = 1.01, fontsize=25)
        
        if i == 6:
            NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
            h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        else:
            h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        
        im = ax.imshow(h_map[i], vmin=h_min, vmax=h_max)
        plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
        plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=17, pad=5) 
        plt.tight_layout()
        #plt.savefig('%s/%s_%i_heatmap_%s.jpeg' % (path_name, date, wavelength, names[i]))
        plt.savefig('%s/%s_%i_heatmap_%s.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
        
        
        flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
    
        fig = plt.figure(figsize=(15,9))
        plt.title('SDO AIA %i.0 Angstrom %s  [Histogram - %s]' % (wavelength, date_title, titles[i]), y = 1.01, fontsize=25)
        plt.xlabel('%s' % cbar_labels[i], fontsize=20, labelpad=10)
        plt.ylabel('Bin Count', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        plt.xlim(h_min, h_max)
        y, x, _ = plt.hist(flat_param, bins=200, range=(h_min, h_max))
        plt.ylim(0, y.max()*1.1)
        #plt.hist(flatten_slopes, bins='auto')  # try this (actually think we want constant bins throughout wavelengths)
        #plt.savefig('%s/%s_%i_Histogram_Slopes.jpeg' % (path_name, date, wavelength))
        plt.savefig('%s/%s_%i_Histogram_%s.pdf' % (path_name, date, wavelength, names[i]), format='pdf')
    
   
    titles_vis = ['Average', 'Middle-File']
    names_vis = ['average', 'mid']
    
    vis = visual
    
    for i in range(2):
        
        fig = plt.figure(figsize=(15,9))
        ax = plt.gca()
        plt.title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles_vis[i]), y = 1.01, fontsize=25)
        #im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
        im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength)
        plt.xlabel('X-position (i) [pixels]', fontsize=20, labelpad=10)
        plt.ylabel('Y-position (j) [pixels]', fontsize=20, labelpad=10)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('Intensity', size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=17, pad=5) 
        plt.tight_layout()
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
        arcsecond coordinate to be converted (int/float)
        
    x2/y1/y2 : 
        Same as x1
    
    image : 
        .fits file to generate subregion image from. (String)
      
    Example:
    ::
        fm.arc2pix(x1,x2,y1,y2, image = 'F:/SDO/data/20130815/171/aia_lev1_171a_2013_08_15t05_59_59_34z_image_lev1.fits.fits') 
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
        pixel coordinate to be converted (int/float)
        
    x2/y1/y2 : 
        Same as x1
    
    image : 
        .fits file to generate subregion image from. (String)
      
    Example:
    ::
        fm.pix2arc(x1,x2,y1,y2, image = 'F:/SDO/data/20130815/171/aia_lev1_171a_2013_08_15t05_59_59_34z_image_lev1.fits.fits') 
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

# add - if rebin = 1, don't add to filename
# have visual / time be both at beginning or end of filename

# add printout of region dimensions and rebinned region dimensions

# update 1/14 - changed 'nf' to % nf - so will print actual # 

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
        fm.datacube(directory='F:/SDO/data/20130530/1600', date='20130530', wavelength=1600,
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

# - insane, I did profiling and 90% of time was spent on finding pixel median values
# - however, it was the median of one pixel.  I took that out and the results are identical
# - in 10% of the time

# t_interp issue -- ran through full program and was exactly same 
# changed t_interp to linspace with one extra point, took out conversion to float in loop

# update 1/14: changed estimated time to reset (not go crazy negative if in same kernel as run before)

# should any of the arrays be reset within loop?

# update 1/16: including 3x3 pixel-box averaging within this function, so fitting function it by itself

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
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate


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
        fm.fft_avg(datacube = DATA, timeseries = TIME, num_seg = 3)
    """
    
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
            x2_box = 2+ii  # if want to use median
            y1_box = 0+jj
            y2_box = 2+jj  # if want to use median
            
            for k in range(0,DATA.shape[0]):
                        im=DATA[k]												# get image
                        #pixmed[k]=np.median(im[x1_box:x2_box,y1_box:y2_box])	# median
                        pixmed[k]= im[x1_box,y1_box]	# median  <-- use this
                        
                        
            # The derotation introduces some bad data towards the end of the sequence. This trims that off
            bad = np.argmax(pixmed <= 0.)		# Look for values <= zero
            last_good_pos = bad - 1			# retain only data before the <=zero
            
            # Get t and v data
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
                 #sig_fft = fftpack.rfft(sig)  # real
                 #sig_fft = np.fft.rfft(sig)  # significantly slower than scipy                 
                 #sig_fft = accelerate.mkl.fftpack.fft(sig)  # possibly use this
                 #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # or this
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
# ran through 1600 rebin4 region and results were identical

# structure for future parallelization - have loop that computes all 3x3 pixel boxes (spectra array variable)
# then run loop through that of curve-fitting, so no dependencies

# update 1/14: changed estimated time to reset (not go crazy negative if in same kernel as run before)

# think will change slope min to 0.3, keep max at 4.0

# think change pl_A from (-0.1, 0.1) to (-0.001 0.001) for M1, and (-0.002, 0.002) for M2
# think change pl_C from (-0.1, 0.1) to (-0.01, 0.01)
# should have changing bounds?  certain wavelengths have higher/lower max/min bounds for certain parameters.

# either change this one, or make another function that does the 3x3 averaging first, then the fits

# removed the 3x3 pixel-box averaging - only fitting now -- pretty sure changed all indices to correct values

# could I keep default method for M1 fit?  

# ran through fitting routine with one region - same results - great

import numpy as np
import scipy.signal
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from matplotlib.widgets import  RectangleSelector
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
import h5py
from scipy import fftpack  # doesnt work in module when called here???
#from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib.pylab as plt
from astropy.convolution import convolve, Box1DKernel
from matplotlib import cm
from numpy.random import randn
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from timeit import default_timer as timer


def spec_fit(spectra_array):
    """
    Calculates pixel-box averages and fits spectra.
    
    spectra_array : 
        3D array of FFT-segment and 3x3-pixel-box averaged region  (array)
        
    Returns : 
        Parameter and M2-fit arrays
      
    Example:
    ::
        fm.spec_fit(spectra_array = SPECTRA)
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
    
    print "The region size is %ii x %ij" % (SPECTRA.shape[0], SPECTRA.shape[1])
    print "%i frequencies were evaluated in the FFT" % SPECTRA.shape[2] 
    
    num_freq = SPECTRA.shape[2]  # determine nubmer of frequencies that are used
        
    # determine frequency values that FFT will evaluate
    freq_size = ((num_freq)*2) + 1  # determined from FFT-averaging script
    time_step = 12  # add as argument, or leave in as constant?
    #time_step = 24  # add as argument, or leave in as constant?
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    
    # initialize arrays to hold parameter values, also each pixel's combined model fit - for tool
    diffM1M2 = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1]))  # dont really use - get rid of?
    params = np.zeros((7, SPECTRA.shape[0], SPECTRA.shape[1]))
    #M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], (len(freqs)+1)/2))  # would save storage / memory space
    M2_fit = np.zeros((SPECTRA.shape[0], SPECTRA.shape[1], SPECTRA.shape[2]))
    
    Uncertainties = np.zeros((6, SPECTRA.shape[0], SPECTRA.shape[1]))
    
    start = timer()
    T1 = 0
    
    
    ### calculate 3x3 pixel-box geometric average.  start at 1 and end 1 before to deal with edges.
    ## 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)
    
    for l in range(0,SPECTRA.shape[0]):
    #for l in range(0,2):
        #print l
        for m in range(0,SPECTRA.shape[1]):
        #for m in range(65,75):
            
                                            
            f = freqs  # frequencies
            s = spectra_array[l][m]  # fourier power
            
            #ds = (1./f**2.2)/1000
            ds = s*0.1  # set the error / variance estimate to a constant percentage of the spectra power-values
            
            # create points to fit model with final parameters 
            #f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],(len(freqs)+1)/2)  # would save storage / memory space
            f_fit = freqs       
        
                    
                    
            ### fit data to models using SciPy's Levenberg-Marquart method
            
            try:
                # initial guesses for fitting parameters
                #P0 = [0.000, 2.0, 0.00003, 0.0022, -6.5, 0.5]
                M1_low = [-0.002, 0.3, -0.01]
                #M1_low = [-0.1, 0.1, -0.1]  # test on 193 - coronal hole
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
            if m > 1:
                P0 = [pl_A[l-1][m-2], slopes[l-1][m-2], pl_C[l-1][m-2], gauss[l-1][m-2], gauss_loc[l-1][m-2], gauss_wid[l-1][m-2]]
                #P0 = [0.0, 1.02, 0.001, 0.001, -4.68, 0.79]
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
            
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
                    
            elif m == 1 and l > 1:
                P0 = [pl_A[l-2][m-1], slopes[l-2][m-1], pl_C[l-2][m-1], gauss[l-2][m-1], gauss_loc[l-2][m-1], gauss_wid[l-2][m-1]]
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
            
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
                    
            else:        
            #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0=[nlfit_l[0], nlfit_l[1], nlfit_l[2], nlfit_g[0], nlfit_g[1], nlfit_g[2]], sigma=ds)
                try:
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)
                    nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=([-0.5,0.1,-0.5,0.00001,-8.5,0.05], [0.5,4.,0.5,0.2,-4.6,0.9]), sigma=ds)
            
                except RuntimeError:
                    print("Error M2 - curve_fit failed")
                
                except ValueError:
                    print("Error M2 - inf/NaN - %i, %i" % (l,m))
            """     
            
            ## fit data to combined power law plus gaussian component model
            #"""        
            try:
                #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(GaussPowerBase, f, s, p0 = P0, bounds=([-4.34,0.5,-8.68,0.00001,-6.5,0.05], [2.,6.,2.,0.2,-4.6,0.8]), sigma=ds)                                  
                M2_low = [-0.002, 0.3, -0.01, 0.00001, -6.5, 0.05]
                #M2_low = [-0.1, 0.1, -0.1, 0.00001, -6.5, 0.05]  # test on 193 - coronal hole
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
            
            
            uncertainties_arr = [dA2, dn2, dC2, dP2, dfp2, dfw2]  # not sure if want to keep these
            Uncertainties[:, l, m] = uncertainties_arr
            
            
            # create model functions from fitted parameters
            m1_fit = PowerLaw(f_fit, A, n, C)
            m2_fit = GaussPowerBase(f_fit, A2,n2,C2,P2,fp2,fw2)
            s_fit_gp_full = GaussPowerBase(f, A2,n2,C2,P2,fp2,fw2)
            m2P_fit = PowerLaw(f_fit, A2, n2, C2)
            m2G_fit = Gauss(f_fit, P2, fp2, fw2)
            
            diffM1M2_temp = (m2_fit - m1_fit)**2  # differences squared
            diffM1M2[l][m] = np.sum(diffM1M2_temp)  # sum of squared differences 
            
            residsgp = (s - s_fit_gp_full)
            redchisqrgp = ((residsgp/ds)**2).sum()/float(f.size-6)
            
            # populate array with parameters
            params[0][l][m] = A2
            params[1][l][m] = n2
            params[2][l][m] = C2
            params[3][l][m] = P2
            params[4][l][m] = fp2
            params[5][l][m] = fw2
            params[6][l][m] = redchisqrgp
            
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
# parameter masking tool
############################
############################
"""     
    
# update 1/13: took out spectra plot on side - not sure if want

# improved formatting a bit - copied over heatmap tool colorbar setup + percentiles

# possibly add a text input for vmin and vmax?

# add the different parameter heatmaps as buttons at the top 
# - will require a bunch of generalization - for the vmin/vmax 

# added to module 1/13, 11:41 PM

# doesn't work - wont stay connected to plot as function call

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from pylab import *
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import figure, show

def mask_param(heatmaps):
    """
    Choose custom vmin or vmax and mask rest of parameter heatmap image.
    
    heatmaps : 
        3D array of parameters  (array)
          
    Example:
    ::
        fm.mask_param(heatmaps = HEATMAPS)
    """
    def update_min(val):
        vmin = s_vmin.val
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i][j] = h_map[1][i][j]
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                if R[i][j] < vmin:
                    R[i][j] = np.nan
     
        ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
        plt.subplots_adjust(bottom=0.25)
        plt.subplots_adjust(left=0.05)
        ax1.set_xlim(0, h_map.shape[2]-1)
        ax1.set_ylim(0, h_map.shape[1]-1)   
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
        #im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # before revision
        im = ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
        cbar.ax.tick_params(labelsize=13, pad=3)
        
    def update_max(val):
        vmax = s_vmax.val
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i][j] = h_map[1][i][j]
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                if R[i][j] > vmax:
                    R[i][j] = np.nan
    
        ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
        plt.subplots_adjust(bottom=0.25)
        plt.subplots_adjust(left=0.05)
        ax1.set_xlim(0, h_map.shape[2]-1)
        ax1.set_ylim(0, h_map.shape[1]-1)  
        ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
        #im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # before revision
        im = ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax)
        cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
        cbar.ax.tick_params(labelsize=13, pad=3)    
        
    
    def reset(event):
        s_vmin.reset()
        s_vmax.reset()
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i][j] = h_map[1][i][j]
        ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)    
        
        
    
    #h_map = np.load('F:/Users/Brendan/Desktop/SolarProject/M2_Spectra_Params/param_20130530_1600_0_296_0_634_numpy.npy')
    h_map = heatmaps
    
    R = np.zeros((h_map.shape[1],h_map.shape[2]))
    
    wavelength = 211
    year = 2012
    month = 9
    day = 23
    date_title = '%i-%02i-%02i' % (year,month,day)
            
    # create list of titles and colorbar names for display on the figures
    titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Power Law Tail', 'Gaussian Amplitude', 'Gaussian Location', 'Gaussian Width', '$/chi^2$', 'Visual Image - Averaged']
    cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location (e^(Value))', 'Width', '$/chi^2$', 'Intensity']
    
    for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i][j] = h_map[1][i][j]
                
    
    # create figure with heatmap and spectra side-by-side subplots
    fig1 = plt.figure(figsize=(20,10))
    ax1 = plt.gca()
    ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1)
    plt.subplots_adjust(bottom=0.25)
    plt.subplots_adjust(left=0.05)
    ax1.set_xlim(0, h_map.shape[2]-1)
    ax1.set_ylim(0, h_map.shape[1]-1)  
    ax1.set_title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[1]), y = 1.01, fontsize=17)
                      
    
    param = h_map[1]  # set initial heatmap to power law index     
    h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    im = ax1.imshow(R, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
    
    #fig, ax = plt.subplots()
    vmin_initial = h_min
    vmax_initial = h_max
    #im = ax1.imshow(R, vmin=0, vmax=2.4, cmap='jet', interpolation='nearest', picker=True)  # gauss amp
    
    # design colorbar for heatmaps
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
    cbar.ax.tick_params(labelsize=13, pad=3)   
            
    
    #ax2 = plt.subplot2grid((1,11),(0, 6), colspan=5, rowspan=1)
    #ax2.loglog()
    #ax2.set_xlim(10**-4.5, 10**-1.3)
    #ax2.set_ylim(10**-5, 10**0)  
    
    axcolor = 'white'
    ax_vmin = plt.axes([0.2, 0.13, 0.45, 0.04], axisbg=axcolor)
    s_vmin = Slider(ax_vmin, 'vmin', 0.1, 3., valinit=vmin_initial)
    
    ax_vmax = plt.axes([0.2, 0.08, 0.45, 0.04], axisbg=axcolor)
    s_vmax = Slider(ax_vmax, 'vmax', 0.1, 3., valinit=vmax_initial)
        #fig.canvas.draw_idle()
    s_vmin.on_changed(update_min)
    s_vmax.on_changed(update_max)
    
    
    resetax = plt.axes([0.72, 0.105, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    
    button.on_clicked(reset)
    
    plt.show()