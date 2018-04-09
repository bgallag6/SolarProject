# -*- coding: utf-8 -*-
"""
Created on Sun Mar 25 21:16:11 2018

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm
from scipy import stats
import matplotlib
import scipy
    
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    #return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2) + C
    return P*np.exp(-0.5*(((np.log(f))-fp)/fw)**2)


#directory = 'F:'
#date = '20130626'
#wavelength = 1600

directory = 'S:'
date = '20111210'
wavelength = 171

#del matplotlib.font_manager.weight_dict['roman']
#matplotlib.font_manager._rebuild()

#matplotlib.rc('text', usetex = True)  # use with latex commands

# 11-param-list
titles = [r'Power Law Slope-Coefficient [flux] - A', r'(b) Power Law Index $n$', r'Power Law Tail - C', r'Gaussian Amplitude [flux] - α', r'(c) Gauss. Loc. $\beta$ [min]', r'(d) Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - α', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
#titles = [r'Power Law Slope-Coefficient [flux] - $A$', r'(b) Power Law Index $n$', r'Power Law Tail - C', r'Gaussian Amplitude [flux] - $\alpha$', r'(c) Gauss. Loc. $\beta$ [min]', r'Gaussian Width - $\sigma$', 'F-Statistic', r'Gaussian Amplitude Scaled - $\alpha$', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
cbar_labels = ['Slope Coefficient', 'Index Value', 'Tail Value', 'Amplitude', 'Location [min]', 'Width', 'F-Statistic', 'Amplitude Scaled', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
names = ['slope_coeff', 'index', 'tail', 'lorentz_amp', 'lorentz_loc', 'lorentz_wid', 'f_test', 'lorentz_amp_scaled', 'r_value', 'roll_freq', 'chisqr']

# load parameter array and visual images from file tree structure 
heatmaps = np.load('%s/DATA/Output/%s/%i/check_if_matches_this/param.npy' % (directory, date, wavelength))
visual = np.load('%s/DATA/Output/%s/%i/visual.npy'% (directory, date, wavelength))  

#visual = visual[1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
visual = visual[1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
h_map = heatmaps    

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels
#plt.rc('font',**{'family':'serif','serif':['Times']})

wavelength = wavelength    

"""
# trim x/y dimensions equally so that resulting region is 1600x1600    
trim_y = int((h_map.shape[1]-1600)/2)
trim_x = int((h_map.shape[2]-1600)/2)
h_map = h_map[:, trim_y:h_map.shape[1]-trim_y, trim_x:h_map.shape[2]-trim_x]  # trim to 1600x1600 (derotate based on mid-file, take off even amounts from both sides)    

x_ticks = [0,200,400,600,800,1000,1200,1400,1600]
y_ticks = [0,200,400,600,800,1000,1200,1400,1600]  
x_ind = [-800,-600,-400,-200,0,200,400,600,800]
y_ind = [800,600,400,200,0,-200,-400,-600,-800]    
"""

xdim = int(np.floor(h_map.shape[2]/100))
ydim = int(np.floor(h_map.shape[1]/100))

x_ticks = [100*i for i in range(xdim+1)]
y_ticks = [100*i for i in range(ydim+1)]

x_ind = x_ticks
y_ind = y_ticks

# generate p-value heatmap + masked Gaussian component heatmaps
df1, df2 = 3, 6  # degrees of freedom for model M1, M2
p_val = ff.sf(h_map[6], df1, df2)


mask_thresh = 0.005  # significance threshold - masked above this value

# mask the Lorentzian component arrays with NaNs if above threshold 
loc_mask = np.copy(h_map[4])
wid_mask = np.copy(h_map[5])    

loc_mask[p_val > mask_thresh] = np.NaN
wid_mask[p_val > mask_thresh] = np.NaN



# determine percentage of region masked 
count = np.count_nonzero(np.isnan(loc_mask))   
total_pix = p_val.shape[0]*p_val.shape[1]
mask_percent = ((np.float(count))/total_pix)*100
            
loc_mask = (1./np.exp(loc_mask))/60.  # convert Gaussian location to minutes


fig_width = 12
fig_height = 10

#for i in range(len(titles)-1):
for i in [1,4,5]:  # use this from now on
#for i in range(4,5):
    
    #fig = plt.figure(figsize=(13,9))
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    #plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
    plt.title('%s' % (titles[i]), y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength
     
    if i == 1:
        h_min = np.percentile(h_map[i],1)  # set heatmap vmin to 1% of data
        h_max = np.percentile(h_map[i],99)  # set heatmap vmax to 99% of data
        cmap = cm.get_cmap('jet', 10)
    elif i == 4:
        h_map[i] = (1./(np.exp(h_map[i])))/60.
        if wavelength in [1600,1700]:
            h_min = 3.
            h_max = 5.
        else:
            h_min = 1.0
            h_max = 11.0
        cmap = cm.get_cmap('jet_r', 10)  # reverse color-scale for Gaussian Location, because of flipped frequencies to seconds
    else:
        flat_param3 = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
        flat_param3 = flat_param3[~np.isnan(flat_param3)]
        #h_min = np.percentile(h_map[i],1)
        #h_max = np.percentile(h_map[i],99)
        h_min = np.percentile(flat_param3,1)
        h_max = np.percentile(flat_param3,99)
        #cmap = 'jet'
        cmap = cm.get_cmap('jet', 10)

    # specify colorbar ticks to be at boundaries of segments
    h_range = np.abs(h_max-h_min)
    h_step = h_range / 10.
    c_ticks = np.zeros((11))
    for h in range(11):
        c_ticks[h] = h_min + h_step*h 
        
    im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=h_min, vmax=h_max)

    plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
    plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
   
    if i != 4:
        cbar = plt.colorbar(im,cax=cax, format='%0.2f')

    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)

    
    
    if i == 4:   
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        
        plt.title(r'(c) Lorentz. Loc. $\beta$ [min]; $p$ $<$ %0.3f \textbar\  $f_{masked}$ = %0.1f\%s' % (mask_thresh, mask_percent, '%'), y = 1.02, fontsize=font_size, fontname="Times New Roman")

        cmap = cm.get_cmap('jet_r', 10)
                   
        im = ax.imshow(np.flipud(loc_mask), cmap = cmap, vmin=h_min, vmax=h_max)

        plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
        plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
        ax.tick_params(axis='both', which='major', pad=10)
        divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
        cax = divider.append_axes("right", size="3%", pad=0.07)

        cbar = plt.colorbar(im,cax=cax, format='%0.1f')

        cbar.ax.tick_params(labelsize=font_size, pad=5) 
        cbar.set_ticks(c_ticks)
        
        #plt.savefig('C:/Users/Brendan/Desktop/%s_%i_%s_mask.pdf' % (date, wavelength, names[i]), format='pdf', bbox_inches='tight')
    
    
    elif i == 5:
        flat_param3 = np.reshape(wid_mask, (wid_mask.shape[0]*wid_mask.shape[1]))
        flat_param3 = flat_param3[~np.isnan(flat_param3)]
        h_min = np.percentile(flat_param3,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(flat_param3,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        
        h_range = np.abs(h_max-h_min)
        h_step = h_range / 10.
        c_ticks = np.zeros((11))
        for h in range(11):
            c_ticks[h] = h_min + h_step*h 
        
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        
        plt.title(r'(d) Lorentz. Wid. $\sigma$ [min]; $p$ $<$ %0.3f \textbar\  $f_{masked}$ = %0.1f\%s' % (mask_thresh, mask_percent, '%'), y = 1.02, fontsize=font_size, fontname="Times New Roman")

        cmap = cm.get_cmap('jet', 10)                    
            
        im = ax.imshow(np.flipud(wid_mask), cmap = cmap, vmin=h_min, vmax=h_max)

        plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
        plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
        ax.tick_params(axis='both', which='major', pad=10)
        divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
        cax = divider.append_axes("right", size="3%", pad=0.07)

        cbar = plt.colorbar(im,cax=cax, format='%0.2f')

        cbar.ax.tick_params(labelsize=font_size, pad=5) 
        cbar.set_ticks(c_ticks)
        
        #plt.savefig('C:/Users/Brendan/Desktop/%s_%i_%s_mask.pdf' % (date, wavelength, names[i]), format='pdf', bbox_inches='tight')
        
        
    
    #"""
    flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
    
    #flat_param = np.reshape(plots[i-2], (plots[i-2].shape[0]*plots[i-2].shape[1]))
    if i == 4:
        flat_param2 = flat_param[~np.isnan(flat_param)]
        sigma2 = np.std(flat_param2)
    
        #h_min2 = np.percentile(flat_param, 1)
        #h_max2 = np.percentile(flat_param, 99)
        
        # calculate some statistics
        mean = np.mean(flat_param2)
        sigma = np.std(flat_param2)   
        
        fig = plt.figure(figsize=(fig_width+1,fig_height))
        plt.title('%i \AA\\ | Lorentzian Location' % wavelength, y = 1.02, fontsize=font_size)  # no date / wavelength
        plt.xlabel('%s' % cbar_labels[i], fontsize=font_size, labelpad=10)
        plt.ylabel('Bin Count', fontsize=font_size, labelpad=10)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        plt.xlim(h_min, h_max)
        y, x, _ = plt.hist(flat_param2, bins=200, range=(h_min, h_max), edgecolor='black')
        
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
                print(gauss_wid)
            
                plt.plot(f,s, linewidth=1.5)
                plt.plot(f,g_fit, linestyle='dashed', linewidth=2.)
                plt.vlines(gauss_center,0,y.max()*1.1, linestyle='dashed', color='red', linewidth=2., label='center=%0.4f' % gauss_center)
            
        
        
        plt.vlines(bin_max, 0, y.max()*1.1, color='black', linestyle='dotted', linewidth=2., label='mode=%0.4f' % bin_max)  
        #plt.vlines(mean, 0, y.max()*1.1, color='red', linestyle='solid', linewidth=1.5, label='mean=%0.6f' % mean)     
        #plt.hlines(y[nearest],fwhm_min,fwhm_min+fwhm, linestyle='dashed', linewidth=2., color='white')
        plt.vlines(0, 0, y.max()*1.1, color='white', linestyle='dashed', linewidth=1.5, label='sigma=%0.4f' % sigma2)
        legend = plt.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
        for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        
        #plt.savefig('C:/Users/Brendan/Desktop/20130626_%i_lorentz_loc_hist.pdf' % wavelength, format='pdf', bbox_inches='tight')
    #"""