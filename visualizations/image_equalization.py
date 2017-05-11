# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 20:04:45 2017

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt    
from skimage import data, img_as_float
from skimage import exposure

directory = 'F:/Users/Brendan/Desktop/SolarProject'
date = '20101208'
wavelength = 1600
#D = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20141227/304/param.npy')
D = np.load('%s/DATA/Output/%s/%i/param.npy' % (directory,date,wavelength))

k = 0
D1 = D[k,50:-50,70:-70]
#D1 = D[k,165:240,350:425] # 20140818
#D1 = D[k,115:200,115:215] # 20140606
#D1 = D[4,90:175,100:210] # 20141227
#D1 = D[4,125:250,250:400] # 20120923
#D1 = D[4,270:375,675:790] # 20130626
#D1 = D[5,250:400,200:300] # 20140818
#D1 = D[1,175:325,260:410] # 20160520
#D1 = D[1,180:420,400:680] # 20131118
#D1 = D[1,40:125,50:145] # 20131118
#D1 = D[k,50:-80,100:-100] # 20121018
#D1 = D[k,135:-165,175:-175] # 20121018

#param1 = H[p1,180:-175,350:-80] # for 20140818
#param2 = H[p2,180:-175,350:-80] # for 20140818
#param1 = H[p1,30:160,40:170] # for 20140819 -- 40:135,55:160
#param2 = H[p2,30:160,40:170] # for 20140819
#param1 = H[p1,70:-55, 190:270] # for 20160426
#param2 = H[p2,70:-55, 190:270] # for 20160426
#param1 = H[p1,160:-130, 260:-165] # for 20160520
#param2 = H[p2,160:-130, 260:-165] # for 20160520
#param1 = H[p1,180:400,390:625] # for 20131118
#param2 = H[p2,180:400,390:625] # for 20131118
#param1 = H[p1,105:205,110:210] # for 20140606
#param2 = H[p2,105:205,110:210] # for 20140606
#param1 = H[p1,60:150,60:150] # for 20140112
#param2 = H[p2,60:150,60:150] # for 20140112
#param1 = H[p1, 325:-250, 350:-125] # for 20160414
#param2 = H[p2, 325:-250, 350:-125] # for 20160414

#param1 = H[p1,150:600,150:650] # for 20170329
#param2 = H[p2,150:600,150:650] # for 20170329
#param1 = H[p1,375:475,625:750] # for 20141025
#param2 = H[p2,375:475,625:750] # for 20141025
#param1 = H[p1,170:320,170:320] # for 20130626
#param2 = H[p2,170:320,170:320] # for 20130626

if k == 4:
    D1 *= -1
#for i in range(21):
#    for j in range(5,6):
        
#        if j == 4:
#            D1 = D[i,j,40:125,50:145]*-1
#        elif j == 2:
#            pass
#        else:
#            #D1 = D[i,j,60:105,70:125] # zoomed
#            D1 = D[i,j,40:125,50:145] # zoomed
        
        #D1 = D1*-1.
"""
        def image_histogram_equalization(image, number_bins=256):
            
            # get image histogram
            image_histogram, bins = np.histogram(image.flatten(), number_bins, normed=True)
            cdf = image_histogram.cumsum() # cumulative distribution function
            cdf = 255 * cdf / cdf[-1] # normalize
        
            # use linear interpolation of cdf to find new pixel values
            image_equalized = np.interp(image.flatten(), bins[:-1], cdf)
            
            #print cdf
            plt.figure()
            plt.hist(cdf,bins=256)
        
            return image_equalized.reshape(image.shape), cdf
        
        
        
        # loop over them
        data_equalized = np.zeros(D1.shape)
        for i in range(1):
            image = D1[:, :]
            data_equalized[:, :] = image_histogram_equalization(image)[0]
            
        plt.figure()
        plt.imshow(D1, vmin=0.4, vmax=1.7)    
        plt.colorbar()
            
        plt.figure()
        plt.imshow(data_equalized)
        plt.colorbar()
"""
        
        
        
        
        
def plot_img_and_hist(img, axes, bins=256):
    """Plot an image along with its histogram and cumulative histogram.

    """
    img = img_as_float(img)
    ax_img, ax_hist = axes
    ax_cdf = ax_hist.twinx()

    # Display image
    ax_img.imshow(img, cmap=plt.cm.gray)
    ax_img.set_axis_off()

    # Display histogram
    ax_hist.hist(img.ravel(), bins=bins, histtype='step', color='black')
    ax_hist.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    ax_hist.set_xlabel('Pixel intensity')
    ax_hist.set_xlim(0, 1)
    ax_hist.set_yticks([])

    # Display cumulative distribution
    img_cdf, bins = exposure.cumulative_distribution(img, bins)
    ax_cdf.plot(bins, img_cdf, 'r')
    ax_cdf.set_yticks([])

    return ax_img, ax_hist, ax_cdf


# Load an example image
img = D1
img /= img.max()

# Contrast stretching
p2 = np.percentile(img, 2)
p98 = np.percentile(img, 98)
img_rescale = exposure.rescale_intensity(img, in_range=(p2, p98))

# Equalization
img_eq = exposure.equalize_hist(img)

# Adaptive Equalization
img_adapteq = exposure.equalize_adapthist(img, clip_limit=0.03)

# Display results
f, axes = plt.subplots(2, 4, figsize=(20, 10))
#plt.suptitle(r'Image / Histogram Equalization Comparison : 20140818 1600$\AA$ Gaussian Amplitude', fontsize=20)
ax_img, ax_hist, ax_cdf = plot_img_and_hist(img, axes[:, 0])
ax_img.set_title('Low contrast image')

y_min, y_max = ax_hist.get_ylim()
ax_hist.set_ylabel('Number of pixels')
ax_hist.set_yticks(np.linspace(0, y_max, 5))

ax_img, ax_hist, ax_cdf = plot_img_and_hist(img_rescale, axes[:, 1])
ax_img.set_title('Contrast stretching')

ax_img, ax_hist, ax_cdf = plot_img_and_hist(img_eq, axes[:, 2])
ax_img.set_title('Histogram equalization')

ax_img, ax_hist, ax_cdf = plot_img_and_hist(img_adapteq, axes[:, 3])
ax_img.set_title('Adaptive equalization')

ax_cdf.set_ylabel('Fraction of total intensity')
ax_cdf.set_yticks(np.linspace(0, 1, 5))

# prevent overlap of y-axis labels
plt.subplots_adjust(wspace=0.4)
plt.show()
#plt.savefig('C:/Users/Brendan/Desktop/1600_segments_20140818/equalizations_20140818_param_%i_seg_%i_zoom.jpeg' % (j,i))
#plt.close()
#plt.savefig('C:/Users/Brendan/Desktop/4_22/image_hist_equalizations_20130626_1600_gauss_loc.pdf', format='pdf')