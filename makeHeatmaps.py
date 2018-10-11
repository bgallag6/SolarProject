# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 11:07:45 2018

@author: Brendan
"""

"""
############################
############################
# generate heatmaps
############################
############################
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm
from scipy import stats
import matplotlib
import yaml

def heatmap(directory, date, wavelength, savefig):
    """
    Generates heatmaps and histograms for each of the parameters:
    power law slope coefficient, power law index, tail value,
    Gaussian amplitude, Gaussian location, Gaussian width, F-statistic,
    p-value masking, Gaussian-component masked plots, rollover period,
    and visual image. 
    
    directory : 
        Directory which contains the DATA folder of the file structure. (string)
        
    date : 
        Date of dataset in 'YYYYMMDD' format. (String)
    
    wavelength :
        Wavelength of dataset. (Integer)
                   
    Example:
    ::
        heatmap(directory = 'S:', date = '20130626', wavelength = '1600', savefig = False) 
    """
    
    # define Gaussian-fitting function
    def Gauss(f, P, fp, fw):
        return P*np.exp(-0.5*((f-fp)/fw)**2)
    
    #matplotlib.rc('text', usetex = True)  # use with latex commands
    
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
    ### for Global EUV Paper
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
    
    """
    xdim = int(np.floor(h_map.shape[2]/100))
    ydim = int(np.floor(h_map.shape[1]/100))
    
    x_ticks = [100*i for i in range(xdim+1)]
    y_ticks = [100*i for i in range(ydim+1)]
    
    x_ind = x_ticks
    y_ind = y_ticks
    """    
    
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
    
    """
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
    """ 
    
    
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
    fig_height = 10  # works better for 20130626
    
    for i in range(h_map.shape[0]):  # use this from now on
    #for i in range(4,5):
        
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
        plt.xticks(x_ticks,x_ind,fontsize=font_size)
        plt.yticks(y_ticks,y_ind,fontsize=font_size)
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
            
            plt.title(r'%s; p < %0.3f | f$_{masked}$ = %0.1f%s' % (titles[i], mask_thresh, mask_percent, '%'), y = 1.02, fontsize=font_size)
            #plt.title(r'(d) %i $\AA$ | f_{masked} = %0.1f%s' % (wavelength, mask_percent, '%'), y = 1.02, fontsize=font_size)
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
            plt.xticks(x_ticks,x_ind,fontsize=font_size)
            plt.yticks(y_ticks,y_ind,fontsize=font_size)
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
    plt.xticks(x_ticks,x_ind,fontsize=font_size)
    plt.yticks(y_ticks,y_ind,fontsize=font_size)
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
    plt.title(r'(d) Rollover Period $T_r$ [min]', y = 1.02, fontsize=font_size)  # no date / wavelength
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
    plt.xticks(x_ticks,x_ind,fontsize=font_size)
    plt.yticks(y_ticks,y_ind,fontsize=font_size)
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
        plt.xticks(x_ticks,x_ind,fontsize=font_size)
        plt.yticks(y_ticks,y_ind,fontsize=font_size)
        ax.tick_params(axis='both', which='major', pad=10)
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="3%", pad=0.07)
        #cbar = plt.colorbar(im,cax=cax)
        #cbar.ax.tick_params(labelsize=font_size, pad=5) 

        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.jpeg' % (directory, date, wavelength, date, wavelength, names_vis[i]))
        #plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.pdf' % (directory, date, wavelength, date, wavelength, names_vis[i]), format='pdf')
        if savefig == True:
            plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_%s.pdf' % (directory, date, wavelength, date, wavelength, names_vis[i]), format='pdf', bbox_inches='tight')
    
        
        # plot sunspot umbra contours on visual image
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
stream = open('specFit_config.yaml', 'r')
cfg = yaml.load(stream)

directory = cfg['temp_dir']
date = cfg['date']
wavelength = cfg['wavelength']
"""

"""
import sys

directory = sys.argv[1]
date = sys.argv[2]
wavelength = int(sys.argv[3])
"""

directory = 'F:'
date = '20130626'
wavelength = 171

heatmap(directory = directory, date= date, wavelength = wavelength, savefig = False)