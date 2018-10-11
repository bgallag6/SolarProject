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

# update 1/21:
# maybe make duplicate fft_avg function, one for computers that have accelerate, one for don't

"""
############################
############################
# generate heatmaps
############################
############################
"""

# think I can take out the trimming of edges?  should be built into the array creation

# when generating heatmaps for rebinned regions - should it scale x/y axes by that factor?

# update 1/19:
# might want to set histogram ranges based on curve_fit parameter bounds?
# want them to be all the same ranges? -- or are we just using to pick out fit problems?

# 3/17:
# changed tick marks from center of segments to boundaries
# think for most part the formatting allows values to be roughly accurate
# changed masking from loop to indexing - speed improvement
# changed gauss location values from seconds to minutes


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm
from scipy import stats
import matplotlib

def heatmap(directory, date, wavelength, savefig):
#def heatmap(heatmaps, visual, date, wavelength, path_name):
    """
    Generates heatmaps and histograms for each of the parameters:
    Power Law Slope Coefficient, Power Law Index, Power Law Tail Value,
    Gaussian Amplitude, Gaussian Location, Gaussian Width, F-Statistic,
    p-value masking, Gaussian-component masked plots, Rollover Period,
    and Visual Images. 
    
    directory : 
        Directory which contains the DATA folder of the file structure. (string)
        
    date : 
        Date of dataset in 'YYYYMMDD' format. (String)
    
    wavelength :
        Wavelength of dataset. (Integer)
                   
    Example:
    ::
        ss.heatmap(heatmap = HEATMAPS, visual = VISUAL, date = '20130815',
           wavelength=211, path_name='C:/Users/Brendan/Desktop/PHYS 326') 
    """
    
    # define Gaussian-fitting function
    #def Gauss(f, P, fp, fw, C):
    def Gauss(f, P, fp, fw):
        #return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2) + C
        return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)
    
    # create arrays to store titles for heatmaps, the names to use when saving the files, and colorbar lables
    #titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'(b) Power Law Index $n$', r'Power Law Tail - $C$', r'Gaussian Amplitude [flux] - $\alpha$', r'(c) Gauss. Loc. $\beta$ [min]', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', 'p-Value']
    #titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'(b) Power Law Index $n$', r'Power Law Tail - $C$', r'Gaussian Amplitude [flux] - $\alpha$', r'(c) Gauss. Loc. $\beta$ [min]', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]'] # 10-param-list
    #titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'Power Law Index $n$', r'Power Law Tail - $C$', r'Gaussian Amplitude [flux] - $\alpha$', r'Gauss. Loc. $\beta$ [min] - ', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', 'p-Value']
    #names = ['slope_coeff', 'index', 'tail', 'gauss_amp', 'gauss_loc', 'gauss_wid', 'f_test', 'gauss_amp_scaled', 'r_value', 'roll_freq']
    #cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [min]', 'Width', 'F-Statistic', 'Amplitude Scaled', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]']
    
    #matplotlib.rc('text', usetex = True)  # use with latex commands
    
    # 11-param-list
    titles = [r'Power Law Slope-Coefficient [flux] - A', r'(b) Power Law Index n', r'Power Law Tail - C', r'Gaussian Amplitude [flux] - α', r'(c) Gauss. Loc. β [min]', r'Gaussian Width - σ', 'F-Statistic', r'Gaussian Amplitude Scaled - α', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
    #titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'(b) Power Law Index $n$', r'Power Law Tail - C', r'Gaussian Amplitude [flux] - $\alpha$', r'(c) Gauss. Loc. $\beta$ [min]', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
    cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [min]', 'Width', 'F-Statistic', 'Amplitude Scaled', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
    names = ['slope_coeff', 'index', 'tail', 'lorentz_amp', 'lorentz_loc', 'lorentz_wid', 'f_test', 'lorentz_amp_scaled', 'r_value', 'roll_freq', 'chisqr']
    
    #vmin = [10**-11, 0.5, 10**-6, 10**-6, -6.5, 0.1, 2.]  # think don't need anymore  (or option to set ranges for specific wavelengths?)
    #vmax = [10**-6, 2.5, 0.003, 10**-2, -4.5, 0.8, 15.]  # think don't need anymore
    
    # load parameter array and visual images from file tree structure 
    heatmaps = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory, date, wavelength))
    visual = np.load('%s/DATA/Output/%s/%i/visual.npy'% (directory, date, wavelength))  
    
    #visual = visual[1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
    visual = visual[:,1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
    h_map = heatmaps    
    
    plt.rcParams["font.family"] = "Times New Roman"
    font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels
    
    wavelength = wavelength    
    
    #"""
    # trim x/y dimensions equally so that resulting region is 1600x1600    
    trim_y = int((h_map.shape[1]-1600)/2)
    trim_x = int((h_map.shape[2]-1600)/2)
    h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    
    
    x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
    y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
    x_ind = [-800,-600,-400,-200,0,200,400,600,800]
    y_ind = [800,600,400,200,0,-200,-400,-600,-800]    
    #"""
    
    #"""
    ### for Global Oscillation Paper
    # trim x/y dimensions equally so that resulting region is 1600x1600    
    #trim_y = (h_map.shape[1]-1200)/2
    #trim_x = (h_map.shape[2]-1200)/2
    #h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    
    
    #x_ticks = [0,200,400,600,800,1000,1200]
    #y_ticks = [0,200,400,600,800,1000,1200]  
    #x_ticks = [10,210,410,610,810,1010,1210]
    #y_ticks = [10,210,410,610,810,1010,1210]  
    #x_ind = [-600,-400,-200,0,200,400,600]
    #y_ind = [600,400,200,0,-200,-400,-600]    
    #"""
    
    #h_map = h_map[:, 0:h_map.shape[1]-50, 0:500]  # for 20130626 blobs      
    #x_ticks = [0,100,200,300,400,500]
    #y_ticks = [0,100,200,300,400]   
    #x_ind = [0,100,200,300,400,500]
    #y_ind = [0,100,200,300,400]   
    #y_ticks = [0,100,200,300] # 20120923
    
    """
    xdim = int(np.floor(h_map.shape[2]/100))
    ydim = int(np.floor(h_map.shape[1]/100))
    
    x_ticks = [100*i for i in range(xdim+1)]
    y_ticks = [100*i for i in range(ydim+1)]
    
    x_ind = x_ticks
    y_ind = y_ticks
    """
    #x_ticks = [0,100,200,300]
    #y_ticks = [0,100,200,300]  
    
    
    # generate p-value heatmap + masked Gaussian component heatmaps
    df1, df2 = 3, 6  # degrees of freedom for model M1, M2
    p_val = ff.sf(h_map[6], df1, df2)
    
    #p_mask = np.copy(p_val)
    
    mask_thresh = 0.005  # significance threshold - masked above this value
       
    p_mask = np.copy(p_val)
    amp_mask = np.copy(h_map[3])
    loc_mask = np.copy(h_map[4])
    wid_mask = np.copy(h_map[5])    
    
    # mask the Gaussian component arrays with NaNs if above threshold 
    p_mask[p_val > mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
    amp_mask[p_val > mask_thresh] = np.NaN
    loc_mask[p_val > mask_thresh] = np.NaN
    wid_mask[p_val > mask_thresh] = np.NaN
    
    #"""
    # for creating sunspot umbra + PPV contour overlays from 1600
    v1600 = np.load('%s/DATA/Output/%s/1600/visual.npy'% (directory, date))  
    v1600 = v1600[:,1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
    p1600 = np.load('%s/DATA/Output/%s/1600/param.npy'% (directory, date)) 
    
    #xdif = v1600.shape[2] - p1600.shape[2]
    #ydif = v1600.shape[1] - p1600.shape[1]
    #v1600 = v1600[:,ydif/2:-ydif/2,xdif/2:-xdif/2]
    #print v1600.shape[1], v1600.shape[2]
    #print p1600.shape[1], p1600.shape[2]
    
    h_map = h_map[:,:p1600.shape[1],:p1600.shape[2]]
    visual = visual[:,:v1600.shape[1],:v1600.shape[2]]
    
    #v_mask = np.copy(visual[0])
    p1600_val = ff.sf(p1600[6], df1, df2)
    #p1600_mask = np.copy(p1600_val)
    v_mask = np.copy(v1600[0])
       
    v_mask[p1600_val < mask_thresh] = 1.  # invert mask, set equal to 1. -- so can make contour
    #""" 
    
    
    # determine percentage of region masked 
    count = np.count_nonzero(np.isnan(p_mask))   
    total_pix = p_val.shape[0]*p_val.shape[1]
    mask_percent = ((np.float(count))/total_pix)*100
                    
    loc_mask = (1./np.exp(loc_mask))/60.  # convert Gaussian location to minutes
    plots = [p_mask, amp_mask, loc_mask, wid_mask]  # make array of masked plots to iterate over
    
    """
    if h_map.shape[2] > h_map.shape[1]:
        aspect_ratio = float(h_map.shape[2]) / float(h_map.shape[1])
        fig_height = 10
        fig_width = 10*aspect_ratio
        
    else:
        aspect_ratio = float(h_map.shape[1]) / float(h_map.shape[2])
        #print aspect_ratio
        #fig_width = 10
        fig_width = 10+2  # works better for 20130626 (with no x/y labels)
        fig_height = 10*aspect_ratio  # works better for 20130626
    """
    
    fig_width = 10+2  # works better for 20130626 (with no x/y labels)
    #fig_width = 10+3  # works better for 20130626 (with x/y labels)
    fig_height = 10  # works better for 20130626
    
    #for i in range(len(titles)-1):
    for i in range(h_map.shape[0]):  # use this from now on
    #for i in range(4,5):
        
        #fig = plt.figure(figsize=(13,9))
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
        plt.title('%s' % (titles[i]), y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength
        
        if i == 6:
            NaN_replace = np.nan_to_num(h_map[i])  # NaN's in chi^2 heatmap were causing issue, replace with 0?
            h_min = np.percentile(NaN_replace,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(NaN_replace,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #cmap = 'jet'  # for continuous colorscale
            cmap = cm.get_cmap('jet', 10)  # specify discrete colorscale with 10 intervals 
        elif i == 4:
            h_map[i] = (1./(np.exp(h_map[i])))/60.
            h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #h_min = np.percentile(h_map[i],0.5)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            #h_max = np.percentile(h_map[i],99.5)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #h_min = 3.5
            #h_max = 4.5
            #cmap = 'jet_r'  # reverse color-scale for Gaussian Location, because of flipped frequencies to seconds
            cmap = cm.get_cmap('jet_r', 10)
        elif i == 9:
            h_map[i] = np.nan_to_num(h_map[i])  # deal with NaN's causing issues
            h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            if wavelength == 1600 or wavelength == 1700:
                h_max = np.percentile(h_map[i],99.9)  # for 1600 sunspot rollover - to show gradient 
            else:
                h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #cmap = 'jet'      
            cmap = cm.get_cmap('jet', 10)               
        else:
            h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
            h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
            #cmap = 'jet'
            cmap = cm.get_cmap('jet', 10)

        # specify colorbar ticks to be at boundaries of segments
        h_range = np.abs(h_max-h_min)
        h_step = h_range / 10.
        c_ticks = np.zeros((11))
        for h in range(11):
            c_ticks[h] = h_min + h_step*h 
            
        im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=h_min, vmax=h_max)
        #im = ax.imshow(h_map[i], cmap = cmap, vmin=h_min, vmax=h_max)
        #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
        #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
        #plt.xticks(x_ticks,fontsize=font_size)
        #plt.yticks(y_ticks,fontsize=font_size)
        plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
        plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
        ax.tick_params(axis='both', which='major', pad=10)
        divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
        cax = divider.append_axes("right", size="3%", pad=0.07)
        if i == 0:
            cbar = plt.colorbar(im,cax=cax, format='%0.2e')
        elif i == 1:
            cbar = plt.colorbar(im,cax=cax, format='%0.2f')
        elif i == 2:
            cbar = plt.colorbar(im,cax=cax, format='%0.4f')
        elif i == 3:
            cbar = plt.colorbar(im,cax=cax, format='%0.4f')
        elif i == 4:
            cbar = plt.colorbar(im,cax=cax, format='%0.1f')
        elif i == 5:
            cbar = plt.colorbar(im,cax=cax, format='%0.2f')
        elif i == 6:
            cbar = plt.colorbar(im,cax=cax, format='%0.1f')
        elif i == 7:
            cbar = plt.colorbar(im,cax=cax, format='%0.2f')
        elif i == 8:
            cbar = plt.colorbar(im,cax=cax, format='%0.2f')
        elif i == 9:
            cbar = plt.colorbar(im,cax=cax, format='%0.2f')
        #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
        cbar.ax.tick_params(labelsize=font_size, pad=5) 
        cbar.set_ticks(c_ticks)
        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.jpeg' % (directory, date, wavelength, date, wavelength, names[i]))
        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf')
        if savefig == True:
            plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf', bbox_inches='tight')
        
        
        if i == 3 or i == 4 or i == 5:   
            fig = plt.figure(figsize=(fig_width,fig_height))
            ax = plt.gca()  # get current axis -- to set colorbar 
            
            plt.title(r'%s; p < %0.3f | f$_{masked}$ = %0.1f%s' % (titles[i], mask_thresh, mask_percent, '%'), y = 1.02, fontsize=font_size, fontname="Times New Roman")
            #plt.title(r'(d) %i $\AA$ | f_{masked} = %0.1f%s' % (wavelength, mask_percent, '%'), y = 1.02, fontsize=font_size, fontname="Times New Roman")
            #plt.title(r'%i$\AA$ Gaussian Location [min]: | f_{masked} = %0.1f' % (wavelength, mask_percent), y = 1.02, fontsize=font_size)
            #plt.title(r'%s: $p$ < %0.3f' % (titles[i], mask_thresh), y = 1.02, fontsize=font_size)
            if i == 4:
                cmap = cm.get_cmap('jet_r', 10)
            else:
                cmap = cm.get_cmap('jet', 10)                    
                
            im = ax.imshow(np.flipud(plots[i-2]), cmap = cmap, vmin=h_min, vmax=h_max)
            #im = ax.imshow(plots[i-2], cmap = cmap, vmin=h_min, vmax=h_max)
            #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
            #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
            #plt.xticks(x_ticks,fontsize=font_size)
            #plt.yticks(y_ticks,fontsize=font_size)
            plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
            plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
            ax.tick_params(axis='both', which='major', pad=10)
            divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
            cax = divider.append_axes("right", size="3%", pad=0.07)
            if i == 3:
                cbar = plt.colorbar(im,cax=cax, format='%0.4f')
            elif i == 4:
                cbar = plt.colorbar(im,cax=cax, format='%0.1f')
            elif i == 5:
                cbar = plt.colorbar(im,cax=cax, format='%0.2f')
            #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
            cbar.ax.tick_params(labelsize=font_size, pad=5) 
            cbar.set_ticks(c_ticks)
            #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s_mask_%i.jpeg' % (directory, date, wavelength, date, wavelength, names[i], (1./mask_thresh)))
            #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s_mask_%i.pdf' % (directory, date, wavelength, date, wavelength, names[i], (1./mask_thresh)), format='pdf')
            if savefig == True:
                plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s_mask_%i.pdf' % (directory, date, wavelength, date, wavelength, names[i], (1./mask_thresh)), format='pdf', bbox_inches='tight')
            
        
        #"""
        flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
        
        # calculate some statistics
        mean = np.mean(flat_param)
        sigma = np.std(flat_param)   
        
        fig = plt.figure(figsize=(fig_width+1,fig_height))
        plt.title('%s' % (titles[i]), y = 1.02, fontsize=font_size)  # no date / wavelength
        plt.xlabel('%s' % cbar_labels[i], fontsize=font_size, labelpad=10)
        plt.ylabel('Bin Count', fontsize=font_size, labelpad=10)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        plt.xlim(h_min, h_max)
        y, x, _ = plt.hist(flat_param, bins=200, range=(h_min, h_max))
        
        #plt.xlim(3.5,5.5)
        #y, x, _ = plt.hist(flat_param, bins=200, range=(3.5,5.5))  # possibly use for 1600/1700 so same range
        
        #n, bins, patches = plt.hist(flat_param, bins=200, range=(h_min, h_max))
        n=y[1:-2]
        bins=x[1:-2]
        elem = np.argmax(n)
        bin_max = bins[elem]
        plt.ylim(0, y.max()*1.1)
        
        
        if i == 4:
            f = x[:-1]
            s = y
            #ds = 1./y
            
            if wavelength == 1600 or wavelength == 1700:
                #nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(Gauss, f, s, method='dogbox', max_nfev=10000)     
                nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(Gauss, f, s)       
                #P, fp, fw, C = nlfit_gp  # unpack fitting parameters
                P, fp, fw = nlfit_gp  # unpack fitting parameters          
                #g_fit = Gauss(f, P,fp,fw, C)  
                g_fit = Gauss(f, P,fp,fw)       
                gauss_center = np.exp(fp)
                gauss_wid = np.exp(fw)
            
                plt.plot(f,s, linewidth=1.5)
                plt.plot(f,g_fit, linestyle='dashed', linewidth=2.)
                plt.vlines(gauss_center,0,y.max()*1.1, linestyle='dashed', color='red', linewidth=2., label='center=%0.4f' % gauss_center)
            
        
        
        plt.vlines(bin_max, 0, y.max()*1.1, color='black', linestyle='dotted', linewidth=2., label='mode=%0.4f' % bin_max)  
        #plt.vlines(mean, 0, y.max()*1.1, color='red', linestyle='solid', linewidth=1.5, label='mean=%0.6f' % mean)     
        #plt.hlines(y[nearest],fwhm_min,fwhm_min+fwhm, linestyle='dashed', linewidth=2., color='white')
        plt.vlines(0, 0, y.max()*1.1, color='white', linestyle='dashed', linewidth=1.5, label='sigma=%0.4f' % sigma)
        legend = plt.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        
        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_Histogram_%s.jpeg' % (directory, date, wavelength, date, wavelength, names[i]))
        if savefig == True:
            plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_Histogram_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf', bbox_inches='tight')
        #"""
    
    ## generate p-value masked plot
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    
    plt.title(r'$p-value$ < %0.3f | f$_{masked}$ = %0.1f%s' % (mask_thresh, mask_percent, '%'), y = 1.02, fontsize=font_size, fontname="Times New Roman")
    cmap = cm.get_cmap('jet', 10)
    
    h_min = 0.0
    h_max = 0.005
    h_range = np.abs(h_max-h_min)
    h_step = h_range / 10.
    c_ticks = np.zeros((11))
    for h in range(11):
        c_ticks[h] = h_min + h_step*h
                  
    im = ax.imshow(np.flipud(plots[0]), cmap = cmap, vmin=h_min, vmax=h_max)
    #im = ax.imshow(plots[0], cmap = cmap, vmin=h_min, vmax=h_max)
    #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
    #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
    #plt.xticks(x_ticks,fontsize=font_size)
    #plt.yticks(y_ticks,fontsize=font_size)
    plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
    plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax, format='%0.4f')
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_%s_mask_%i.jpeg' % (directory, date, wavelength, date, wavelength, names[i], (1./mask_thresh)))
    if savefig == True:
        plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_p_value_mask_%i.pdf' % (directory, date, wavelength, date, wavelength, (1./mask_thresh)), format='pdf', bbox_inches='tight')
            
            
    # generate 'rollover frequency' heatmap
    roll_freq = (h_map[2] / h_map[0])**(-1./ h_map[1])
    roll_freq = (1./roll_freq)/60.
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    plt.title(r'(d) Rollover Period $T_r$ [min]', y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength
    roll_freq = np.nan_to_num(roll_freq)  # deal with NaN's causing issues
    h_min = np.percentile(roll_freq,1.)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    if wavelength == 1600 or wavelength == 1700:
        #h_max = np.percentile(roll_freq,99.9)  # for 1600 sunspot rollover - to show gradient 
        h_max = np.percentile(roll_freq,99.99)  # for 1600 sunspot rollover - to show gradient 
    else:
        h_max = np.percentile(roll_freq,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    
    h_range = np.abs(h_max-h_min)
    h_step = h_range / 10.
    c_ticks = np.zeros((11))
    for h in range(11):
        c_ticks[h] = h_min + h_step*h
     
    cmap = cm.get_cmap('jet', 10)    
    
    im = ax.imshow(np.flipud(roll_freq), cmap = cmap, vmin=h_min, vmax=h_max)
    #im = ax.imshow(roll_freq, cmap = cmap, vmin=h_min, vmax=h_max)
    #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
    #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
    #plt.xticks(x_ticks,fontsize=font_size)
    #plt.yticks(y_ticks,fontsize=font_size)
    plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
    plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax,format='%0.1f')
    #cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_roll_freq.jpeg' % (directory, date, wavelength, date, wavelength))
    #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_roll_freqB.pdf' % (directory, date, wavelength, date, wavelength), format='pdf')
    if savefig == True:
        plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_roll_freqB.pdf' % (directory, date, wavelength, date, wavelength), format='pdf', bbox_inches='tight')
    
    
    
    #"""
    # generate visual images
    titles_vis = ['Average', 'Middle-File']
    names_vis = ['average', 'mid']
    
    vis = visual
    #trim_yv = (vis.shape[0]-1600)/2
    #trim_xv = (vis.shape[1]-1600)/2
    #vis = vis[trim_yv:vis.shape[0]-trim_yv, trim_xv:vis.shape[1]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides) 
    
    #trim_yv = (vis.shape[1]-1600)/2
    #trim_xv = (vis.shape[2]-1600)/2
    #vis = vis[:, trim_yv:vis.shape[1]-trim_yv, trim_xv:vis.shape[2]-trim_xv]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    
           
    
    #for i in range(2):
    for i in range(1):
        
        v_min = np.percentile(vis[i],1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        v_max = np.percentile(vis[i],99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)  
        
        #v_min = np.percentile(vis,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        #v_max = np.percentile(vis,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)  
        
        #fig = plt.figure(figsize=(12,9))
        fig = plt.figure(figsize=(fig_width,fig_height))
        
        ax = plt.gca()
        #ax = plt.subplot2grid((1,31),(0, 0), colspan=30, rowspan=1)  #to substitute for colorbar space
        #plt.subplots_adjust(right=0.875)  #to substitute for colorbar space
        plt.title('(a) Visual %s' % (titles_vis[i]), y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength
        #plt.title('(f) %i $\AA$' % wavelength, y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength
        #im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
        #im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)   
        im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
        #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
        #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
        #plt.xticks(x_ticks,fontsize=font_size)
        #plt.yticks(y_ticks,fontsize=font_size)
        plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
        plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
        ax.tick_params(axis='both', which='major', pad=10)
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="3%", pad=0.07)
        #cbar = plt.colorbar(im,cax=cax)
        #cbar.ax.tick_params(labelsize=font_size, pad=5) 

        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.pdf' % (directory, date, wavelength, date, wavelength, names_vis[i]), format='pdf')
        if savefig == True:
            plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.pdf' % (directory, date, wavelength, date, wavelength, names_vis[i]), format='pdf', bbox_inches='tight')
    
    
        dates = ['20101208','20111210','20121018', '20131118', '20140112', '20140606', '20140818', '20140910', '20141227', '20150104','20160414', '20160426', '20160520', '20160905', '20170329']
        contours = [85., 125., 95., 110., 67., 95., 75., 115., 75., 113., 87., 60., 60., 95., 83.]
        
        if date in dates:
            
            # create inverted p-value mask - to plot as contour 
            delta = 1.
            x = np.arange(0., v_mask.shape[1], delta)
            y = np.arange(0., v_mask.shape[0], delta)
            X, Y = np.meshgrid(x, y)
            Z1 = np.flipud(v_mask)    
            
            k = dates.index(date)
            #print contours[k]
    
            #fig = plt.figure(figsize=(12,9))
            fig = plt.figure(figsize=(fig_width,fig_height))
            
            ax = plt.gca()
            plt.title('Visual: %s' % (titles_vis[i]), y = 1.02, fontsize=font_size)  # no date / wavelength
            #im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])  
            im = ax.imshow(np.flipud(vis[i]), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
            CS = plt.contour(X, Y, Z1, levels=[2.], linewidths=2, colors='red', linestyles='solid')
            CS = plt.contour(X, Y, Z1, levels=[2.], linewidths=2, colors='white', linestyles='dotted')
            CS = plt.contour(X, Y, np.flipud(v1600[i]), levels=[contours[k]], linewidths=2, colors='black', linestyles='solid')
            CS = plt.contour(X, Y, np.flipud(v1600[i]), levels=[contours[k]], linewidths=2, colors='white', linestyles='dashed')
            #im = ax.imshow(np.flipud(v_mask), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
            #plt.xlabel('X-Position [Pixels]', fontsize=font_size, labelpad=10)
            #plt.ylabel('Y-Position [Pixels]', fontsize=font_size, labelpad=10)
            #plt.xticks(x_ticks,fontsize=font_size)
            #plt.yticks(y_ticks,fontsize=font_size)
            #plt.xticks(x_ticks,x_ind,fontsize=font_size)
            #plt.yticks(y_ticks,y_ind,fontsize=font_size)
            ax.tick_params(axis='both', which='major', pad=10)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.07)
            cbar = plt.colorbar(im,cax=cax)
            cbar.ax.tick_params(labelsize=font_size, pad=5) 
            #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s_umbra.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
            if savefig == True:
                plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s_umbra.pdf' % (directory, date, wavelength, date, wavelength, names_vis[i]), format='pdf')
    #"""



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
    print("[Arcsecond Coordinate] = [Pixel Coordinate]")
    print("%i arcsec = x1 = %i pixel" % (x1,x1_coord))
    print("%i arcsec = x2 = %i pixel" % (x2,x2_coord))
    print("%i arcsec = y1 = %i pixel" % (y1,y1_coord))
    print("%i arcsec = y2 = %i pixel" % (y2,y2_coord))


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
    print("[Pixel Coordinate] = [Arcsecond Coordinate]")
    print("%i pixel = x1 = %i arcsec" % (x1,x1_coord))
    print("%i pixel = x2 = %i arcsec" % (x2,x2_coord))
    print("%i pixel = y1 = %i arcsec" % (y1,y1_coord))
    print("%i pixel = y2 = %i arcsec" % (y2,y2_coord))    


   
   
"""
#######################
#######################
# datacube + derotate #
#######################
#######################
"""


from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
import numpy as np
import astropy.units as u
from astropy.time import Time



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
        return eval(''.join(evList))
        
    # define subregion coordinates   
    x1 = sub_reg_coords[0]  # j1
    x2 = sub_reg_coords[1]  # j2
    y1 = sub_reg_coords[2]  # i1
    y2 = sub_reg_coords[3]  # i2
    
    # create a list of all the files. This is USER-DEFINED
    flist = sorted(glob.glob('%s/FITS/%s/%i/aia*.fits' % (directory,date,wavelength)))
    #flist = sorted(glob.glob('S:/FITS/%s/%i/aia*.fits' % (date,wavelength)))
    nf = len(flist)

    # Select the image that is the "middle" of our selection.
    # We do this because the solar derotation algorithm operates centered on the 
    mid_file = np.int(np.floor(nf / 2))
    
    mc_list = []  # create an empty list
    
    count = 0  # counter to see where program is at 
    
    # Use defined coordinates, extract the submaps from each AIA image, and store
    # them in the empty list. This takes many minutes to complete.
    print(" ")
    print("Reading files and extracting submaps. This takes a while...")
    print(" ")
   
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
    print(" ")
    print("Creating derotated cube...")
    print("Please wait...")
    
    # Perform the derotation of the submaps. This take a while too.
    #dr = mapcube_solar_derotate(new_mapcube)
    dr = mapcube_solar_derotate(new_mapcube, layer_index=mid_file)
    
    print("done derotating")
    
    
    mid_subarr = dr[mid_file].data		# extract data from middle file of derotated datacube
    rem_i = mid_subarr.shape[0] % bin_frac  # calculate remainder to get integer dimensions
    rem_j = mid_subarr.shape[1] % bin_frac  # calculate remainder to get integer dimensions
    subarr_idim = (mid_subarr.shape[0] - rem_i) / bin_frac  # get ydim
    subarr_jdim = (mid_subarr.shape[1] - rem_j) / bin_frac  # get xdim
     
    mc_list = np.zeros((1,2,3))  # free up memory
    new_mapcube = np.zeros((1,2,3))  # free up memory
    
    #t = dr[0].date  # extract the date / time from the first image
    #base_time = (t.hour * 3600.) + (t.minute * 60.) + t.second  # convert date / time to seconds

    # initialize arrays to hold exposure time, pixel data, and time values
    I = np.empty((nf))  # exposure time
    DATA = np.empty((nf, subarr_idim, subarr_jdim), dtype=np.int16)  # save as int16, since that is what original is
    TIME = np.empty((nf))  # time stamps
    AVG = np.zeros((subarr_idim, subarr_jdim))
    
    # loop through datacube and extract pixel data and time values
    """# this is probably another loop that I could take out and extract directly"""
    for p in range(nf):
        Ex = dr[p].exposure_time
        I[p] = Ex.value
        L = dr[p].data
        L_trim = L[0:(mid_subarr.shape[0] - rem_i), 0:(mid_subarr.shape[1] - rem_j)]
        small_L = rebin(L_trim, L_trim.shape[0]/bin_frac, L_trim.shape[1]/bin_frac)
        DATA[p][:][:] = small_L
        AVG += (small_L / Ex.value)  # create normalized average visual image
        T = dr[p].date
        #curr_time=(T.hour * 3600.)+(T.minute * 60.)+T.second	
        #TIME[p] = curr_time - base_time  # calculate running time of image
        TIME[p] = Time(T).jd  # extract julian day time from each image
        
    #TIME[TIME < 0] += 86400  # for times that cross from 23:59 - 00:01 --> add 1 day to correct for discrepancy
    TIME -= TIME[0]  # set all values equal to time since first entry
    TIME = np.around(TIME*86400)  # get the time value in seconds, and round to nearest whole number
    
    # save the pixel-value, time-array, and exposure-time datacubes as numpy files
    #np.save('%s/DATA/Temp/%s/%i/%i_%ii_%i_%ij_data_rebin%i.npy' % (directory, date, wavelength, y1, y2, x1, x2, bin_frac), DATA)
    #np.save('%s/DATA/Temp/%s/%i/%i_%ii_%i_%ij_time.npy' % (directory, date, wavelength, y1, y2, x1, x2), TIME)
    #np.save('%s/DATA/Temp/%s/%i/%i_%ii_%i_%ij_exposure.npy' % (directory, date, wavelength, y1, y2, x1, x2), I)
    np.save('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength), DATA)
    np.save('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength), TIME)
    np.save('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength), I)
    
    # calculate the average-intensity image of the timeseries 
    #AVG = np.average(DATA,axis=0)
    AVG /= nf
    
    # determine the middle file of the timeseries
    mid_num = (DATA.shape[0]/2)
    mid = DATA[mid_num]
    
    print("Middle file is number %i" % mid_num)
    
    # check the two image sizes agree
    print(" Average Image Dimensions = %i, %i" % (AVG.shape[0], AVG.shape[1]))
    print(" Middle Image Dimensions = %i, %i" % (mid.shape[0], mid.shape[1]))
    
    # store average and middle images in array
    visual = np.zeros((2,AVG.shape[0],AVG.shape[1]))
    visual[0] = AVG
    visual[1] = mid
    
    print(visual.shape)  # check array size agrees with expected
    
    # save visual-image array
    #np.save('%s/DATA/Output/%s/%i/%i_%ii_%i_%ij_visual.npy' % (directory, date, wavelength, y1, y2, x1, x2), visual)
    np.save('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength), visual)
    
    # generate images of each visual region, to see if as expected
    #fig = plt.figure(figsize=(20,20))
    #plt.imshow(visual[0])
    #fig = plt.figure(figsize=(20,20))
    #plt.imshow(visual[1])
    
    
    

"""
###################################################
###################################################
# FFT segment averaging + 3x3 Pixel Box Averaging #
###################################################
###################################################
"""

# 1/29: saving datacube as int, normalizing by exposure time in loop (saves 4x memory)

# 2/23: changed from geometric 3x3 average to arithmetic average


import numpy as np
import scipy.signal
from pylab import *
from sunpy.map import Map
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate  # put inside function


def fft_avg(directory, date, wavelength, num_seg):
    """
    Calculates segment-averaged FFT for region, and 3x3 pixel-box average.
    
    directory : 
        The directory containing the derotated datacube, as well as the time and exposure arrays. (string)
        
    date :
        The date of the timeseries.  (string)
                
    wavelength :
        The wavelength of the dataset.  (int)
        
    num_seg :
        The number of segments to divide the timeseries into.  (int)
        
    Returns : 
        3D array of FFT-segment and 3x3-pixel-box averaged region.
      
    Example:
    ::
        ss.fft_avg(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength, num_seg = 6)
    """
    #import accelerate  # switch on if computer has installed
    
    from scipy import fftpack
    
    DATA = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))
    
    TIME = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))
    
    Ex = np.load('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength))
    
    print(DATA.shape) 
    
    print("Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0]))
          
    ## determine frequency values that FFT will evaluate
    if wavelength == 1600 or wavelength == 1700:
      time_step = 24  # 24-second cadence for these wavelengths
    else:
      time_step = 12  # 12-second cadence for the others
      #time_step = 24  # for half-cadence try
      
    t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/time_step)+1)  # interpolate onto default-cadence time-grid
      
    n_segments = num_seg  # break data into 12 segments of equal length
    n = len(t_interp)
    rem = n % n_segments
    freq_size = (n - rem) // n_segments 
    
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    reslt = (DATA.shape[0] == TIME.shape[0])
    print("DATA and TIME array sizes match: %s" % reslt)
    
    pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
    spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))
    
    print("length time-interp array = %i" % n)
    print("size for FFT to consider = %i" % freq_size)
    print("length of sample freq array = %i" % len(sample_freq))
    print("length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs))
    
    
    start = timer()
    T1 = 0
    
    for ii in range(spectra_seg.shape[0]):
    #for ii in range(0,5):
    
        for jj in range(spectra_seg.shape[1]):
        #for jj in range(0,5):        

            pixmed = DATA[:,ii,jj] / Ex  # extract timeseries + normalize by exposure time   
        
            v_interp = np.interp(t_interp,TIME,pixmed)  # interpolate pixel-intensity values onto specified time grid
            
            data = v_interp
            
            avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
    
            data = data[0:len(data)-rem]  # trim timeseries to be integer multiple of n_segments
            split = np.split(data, n_segments)  # create split array for each segment
       
            for i in range(0,n_segments):               
                
              ## perform Fast Fourier Transform on each segment       
              sig = split[i]
              sig_fft = fftpack.fft(sig)
              #sig_fft = fftpack.rfft(sig)  # real-FFT                
              #sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
              #sig_fft = accelerate.mkl.fftpack.rfft(sig)  # this is slightly faster
              powers = np.abs(sig_fft)[pidxs]
              norm = len(sig)
              powers = ((powers/norm)**2)*(1./(sig.std()**2))*2   # normalize the power
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
            print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spectra_seg.shape[0], T_hr, T_min, T_sec))
        else:
            T_est2 = T2*(spectra_seg.shape[0]-ii)
            T_min2, T_sec2 = divmod(T_est2, 60)
            T_hr2, T_min2 = divmod(T_min2, 60)
            print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spectra_seg.shape[0], T_hr2, T_min2, T_sec2))
        T1 = T
        
        
    # print estimated and total program time to screen 
    print("Beginning Estimated time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec))
    T_act = timer() - start
    T_min3, T_sec3 = divmod(T_act, 60)
    T_hr3, T_min3 = divmod(T_min3, 60)
    print("Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3)) 
    
    
    # initialize arrays to hold temporary results for calculating arithmetic average (changed from geometric)
    temp = np.zeros((9,spectra_seg.shape[2]))  # maybe have 3x3 to be generalized   
    #temp = np.zeros((25,spectra_seg.shape[2]))
    spectra_array = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))
    spectra_StdDev = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))
    #spectra_array = np.zeros((spectra_seg.shape[0]-4, spectra_seg.shape[1]-4, spectra_seg.shape[2]))    
    
    ### calculate 3x3 pixel-box arithmetic average.  start at 1 and end 1 before to deal with edges.
    ## previously for geometric -- 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)

    for l in range(1,spectra_seg.shape[0]-1):
    #for l in range(2,spectra_seg.shape[0]-2):
    #for l in range(1,25):
        #print l
        for m in range(1,spectra_seg.shape[1]-1):
        #for m in range(2,spectra_seg.shape[1]-2):
        #for m in range(1,25):
        
            """
            temp[0] = np.log10(spectra_seg[l-1][m-1])  # previous geometric calculations
            temp[1] = np.log10(spectra_seg[l-1][m])
            temp[2] = np.log10(spectra_seg[l-1][m+1])
            temp[3] = np.log10(spectra_seg[l][m-1])
            temp[4] = np.log10(spectra_seg[l][m])
            temp[5] = np.log10(spectra_seg[l][m+1])
            temp[6] = np.log10(spectra_seg[l+1][m-1])
            temp[7] = np.log10(spectra_seg[l+1][m])
            temp[8] = np.log10(spectra_seg[l+1][m+1])
            """
            
            #"""
            temp[0] = spectra_seg[l-1][m-1]
            temp[1] = spectra_seg[l-1][m]
            temp[2] = spectra_seg[l-1][m+1]
            temp[3] = spectra_seg[l][m-1]
            temp[4] = spectra_seg[l][m]
            temp[5] = spectra_seg[l][m+1]
            temp[6] = spectra_seg[l+1][m-1]
            temp[7] = spectra_seg[l+1][m]
            temp[8] = spectra_seg[l+1][m+1]
            #"""
            
            """            
            w = 0.08
            wc = 1. - (8*w)
            
            temp[0] = spectra_seg[l-1][m-1]*w
            temp[1] = spectra_seg[l-1][m]*w
            temp[2] = spectra_seg[l-1][m+1]*w
            temp[3] = spectra_seg[l][m-1]*w
            temp[4] = spectra_seg[l][m]*wc
            temp[5] = spectra_seg[l][m+1]*w
            temp[6] = spectra_seg[l+1][m-1]*w
            temp[7] = spectra_seg[l+1][m]*w
            temp[8] = spectra_seg[l+1][m+1]*w
            """            

            p_avg = np.average(temp, axis=0)

            spectra_array[l-1][m-1] = p_avg
            spectra_StdDev[l-1][m-1] = np.std(temp, axis=0)
            
            
    np.save('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength), spectra_array)
    np.save('%s/DATA/Temp/%s/%i/3x3_stddev.npy' % (directory, date, wavelength), spectra_StdDev)
    #np.save('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength), spectra_seg)  # for no 3x3 pixel box averaging
    #return spectra_array
    
    
    
    
"""
##########################
##########################
# FFT Overlapping Window #
##########################
##########################
"""


import numpy as np
import scipy.signal
from pylab import *
from sunpy.map import Map
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
#from scipy import fftpack  # not working with this called here???
from timeit import default_timer as timer
#import accelerate  # put inside function


def fft_overlap(directory, date, wavelength, window_length, overlap_pct, pixel_box):
    """
    Calculates segment-averaged FFT for region, and 3x3 pixel-box average.
    
    directory : 
        The directory containing the derotated datacube, as well as the time and exposure arrays. (string)
        
    date :
        The date of the timeseries.  (string)
                
    wavelength :
        The wavelength of the dataset.  (int)
        
    num_seg :
        The number of segments to divide the timeseries into.  (int)
        
    Returns : 
        3D array of FFT-segment and 3x3-pixel-box averaged region.
      
    Example:
    ::
        ss.fft_avg(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength, num_seg = 6)
    """
    import accelerate  # switch on if computer has installed
    
    from scipy import fftpack
    
    DATA = np.load('%s/DATA/Temp/%s/%i/derotated.npy' % (directory, date, wavelength))
    
    TIME = np.load('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength))
    
    Ex = np.load('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength))
    
    print(DATA.shape)
    
    print("Number of seconds in timeseries = %i" % (TIME[len(TIME)-1] - TIME[0]))    
        
    ## determine frequency values that FFT will evaluate
    if wavelength == 1600 or wavelength == 1700:
      time_step = 24  # 24-second cadence for these wavelengths
    else:
      time_step = 12  # 12-second cadence for the others
      
    t_interp = np.linspace(0, TIME[len(TIME)-1], (TIME[len(TIME)-1]/time_step)+1)  # interpolate onto default-cadence time-grid
        
    freq_size = (int(window_length[0:2])*3600 + int(window_length[3:5])*60)/time_step
    #rem = len(t_interp) % 
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)
    freqs = sample_freq[pidxs]
    
    reslt = (DATA.shape[0] == TIME.shape[0])
    print("DATA and TIME array sizes match: %s" % reslt)
    
    #pixmed=np.empty(DATA.shape[0])  # Initialize array to hold median pixel values
    #spectra_seg = np.zeros((DATA.shape[1],DATA.shape[2],len(freqs)))
    
    n_segments = ((len(t_interp)-freq_size) / int(freq_size*(1.-(overlap_pct/100.)))) + 1   
    print(n_segments)
    print(len(freqs))
    spec_array = np.zeros((n_segments,DATA.shape[1],DATA.shape[2],len(freqs)))
    
    print("length time-interp array = %i" % len(t_interp))
    print("size for FFT to consider = %i" % freq_size)
    print("length of sample freq array = %i" % len(sample_freq))
    print("length of freqs array = %i (should be 1/2 of two above rows)" % len(freqs))
    
    
    start = timer()
    T1 = 0
    
    for ii in range(spec_array.shape[1]):
    #for ii in range(0,5):
    
        for jj in range(spec_array.shape[2]):
        #for jj in range(0,5):        
            
            pixmed = DATA[:,ii,jj] / Ex  # extract timeseries + normalize by exposure time     
        
            v_interp = np.interp(t_interp,TIME,pixmed)  # interpolate pixel-intensity values onto specified time grid
            
            avg_array = np.zeros((len(freqs)))  # initialize array to hold fourier powers
            avg_smooth = np.zeros((len(freqs)))  # initialize array to hold fourier powers
            
            sig = np.zeros((len(freqs)))
            #sig_test1 = sig[595:600]
            points_overlap = len(sig)*overlap_pct/100
            points_nonoverlap = (len(sig)-points_overlap)
            #print points_nonoverlap
            for i in range(n_segments):  
                sig = v_interp[i*points_nonoverlap:(i*points_nonoverlap)+freq_size] # i*#overlap points + freqs_size
                #print (i*points_nonoverlap), (i*points_nonoverlap)+600
                #sig_fft = fftpack.fft(sig)
                sig_fft = accelerate.mkl.fftpack.fft(sig)  # MKL-accelerated is (2x) faster
                powers = np.abs(sig_fft)[pidxs]
                norm = len(sig)
                powers = ((powers/norm)**2)*(1./(sig.std()**2))*2  # normalize the power
                avg_array += powers
                      
                orig_powerspec = powers
                # Now Smooth the Power Spectra
                orig = orig_powerspec.copy()		# make a copy            
                window_len=7   # Keep this an odd number            
                window='hanning'  # OTHER CHOICES =  'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                s=np.r_[orig[window_len-1:0:-1],orig,orig[-1:-window_len:-1]]            
                if window == 'flat': #moving average
                	w=np.ones(window_len,'d')
                else:
                	w=eval('np.'+window+'(window_len)')                        	
                y=np.convolve(w/w.sum(),s,mode='valid')
                y=y[( window_len/2+1): len(orig)+( window_len/2 +1)]	# Crop down the new TS            
                # Now smooth the smooth...
                y2c = y.copy()
                s2=np.r_[y2c[window_len-1:0:-1],y2c,y2c[-1:-window_len:-1]]
                y2=np.convolve(w/w.sum(),s2,mode='valid')
                y2=y2[( window_len/2+1): len(orig)+( window_len/2 +1)]   # Crop down the new TS
                avg_smooth += y2
                
                avg_raw = np.copy(avg_array)  # take the average of the segments
                avg_raw = avg_raw / (i+1)
                
                avg_smooth2 = np.copy(avg_smooth)
                avg_smooth2 = avg_smooth2 / (i+1)
                
                #spec_array[0][i] = powers
                spec_array[i][ii][jj] = powers
                #spec_array[i][ii][jj] = y2            
        
        
        # estimate time remaining and print to screen
        T = timer()
        T2 = T - T1
        if ii == 0:
            T_init = T - start
            T_est = T_init*(spec_array.shape[1])  
            T_min, T_sec = divmod(T_est, 60)
            T_hr, T_min = divmod(T_min, 60)
            print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spec_array.shape[1], T_hr, T_min, T_sec))
        else:
            T_est2 = T2*(spec_array.shape[1]-ii)
            T_min2, T_sec2 = divmod(T_est2, 60)
            T_hr2, T_min2 = divmod(T_min2, 60)
            print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (ii, spec_array.shape[1], T_hr2, T_min2, T_sec2))
        T1 = T
        
        
    # print estimated and total program time to screen 
    print("Beginning Estimated time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec))
    T_act = timer() - start
    T_min3, T_sec3 = divmod(T_act, 60)
    T_hr3, T_min3 = divmod(T_min3, 60)
    print("Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3))
    
    if pixel_box == True:
        # initialize arrays to hold temporary results for calculating arithmetic average (changed from geometric)
        temp = np.zeros((9,spec_array.shape[3]))  # maybe have 3x3 to be generalized   
        p_avg = np.zeros((spec_array.shape[3]))  # would pre-allocating help? (seems to)
        spectra_array = np.zeros((spec_array.shape[0], spec_array.shape[1]-2, spec_array.shape[2]-2, spec_array.shape[3]))  # would pre-allocating help? (seems to)
                  
        ### calculate 3x3 pixel-box arithmetic average.  start at 1 and end 1 before to deal with edges.
        ## previously for geometric -- 10^[(log(a) + log(b) + log(c) + ...) / 9] = [a*b*c*...]^(1/9)
        for k in range(spec_array.shape[0]):
            for l in range(1,spec_array.shape[1]-1):
            #for l in range(1,25):
                #print l
                for m in range(1,spec_array.shape[2]-1):
                #for m in range(1,25):
                               
                    temp[0] = spec_array[k][l-1][m-1]
                    temp[1] = spec_array[k][l-1][m]
                    temp[2] = spec_array[k][l-1][m+1]
                    temp[3] = spec_array[k][l][m-1]
                    temp[4] = spec_array[k][l][m]
                    temp[5] = spec_array[k][l][m+1]
                    temp[6] = spec_array[k][l+1][m-1]
                    temp[7] = spec_array[k][l+1][m]
                    temp[8] = spec_array[k][l+1][m+1]
        
        
                    temp9 = np.sum(temp, axis=0)
                    p_avg = temp9 / 9.
                    #spectra_array[l-1][m-1] = np.power(10,p_geometric)
                    spectra_array[k][l-1][m-1] = p_avg
                    
        spec_array = spectra_array
        np.save('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength), spectra_array)
    else:
        np.save('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength), spec_array)
    #return spectra_array
    
    


"""
############################
############################
# Save Memory-Mapped Array
############################
############################
"""
import numpy as np

def mem_map(directory, date, wavelength):
    """
    Creates memory-mapped array of spectra array to use in MPI
    
    directory : 
        The directory containing the spectra array. (string)
        
    date :
        The date of the timeseries.  (string)
                
    wavelength :
        The wavelength of the dataset.  (int)
        
    Returns : 
        Memory-mapped spectra array. 
      
    Example:
    ::
        ss.mem_map(directory='%s' % (directory), date='%s' % (date), wavelength= wavelength)
    """
    
    # load original array 
    original = np.load('%s/DATA/Temp/%s/%i/spectra.npy' % (directory, date, wavelength))
    original_StdDev = np.load('%s/DATA/Temp/%s/%i/3x3_stddev.npy' % (directory, date, wavelength))
    print(original.shape)
    
    if original.ndim == 3:
        orig_shape = np.array([original.shape[0], original.shape[1], original.shape[2]])
        
        # create memory-mapped array with similar datatype and shape to original array
        mmap_arr = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='w+', shape=tuple(orig_shape))
        mmap_StdDev = np.memmap('%s/DATA/Temp/%s/%i/3x3_stddev_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='w+', shape=tuple(orig_shape))
        
        
    elif original.ndim == 4:
        orig_shape = np.array([original.shape[0], original.shape[1], original.shape[2], original.shape[3]])
        
        # create memory-mapped array with similar datatype and shape to original array
        mmap_arr = np.memmap('%s/DATA/Temp/%s/%i/spectra_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='w+', shape=(original.shape[0],original.shape[1],original.shape[2], original.shape[3]))
     
    # write data to memory-mapped array
    mmap_arr[:] = original[:]
    mmap_StdDev[:] = original_StdDev[:]
    
    # flush memory changes to disk, then remove memory-mapped object
    del mmap_arr
    del mmap_StdDev
    
    np.save('%s/DATA/Temp/%s/%i/spectra_mmap_shape.npy' % (directory, date, wavelength), orig_shape)