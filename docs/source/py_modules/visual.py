import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

def visual(dataset, date, wavelength, path_name):
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
        fm.get_data(wavelength=1600, time_begin='2016/09/23 00:00:00', time_end='2016/09/23 00:05:00', path_name='C:/Users/Brendan/Desktop/SDO_test')
    """

    titles = ['Average', 'Middle-File']
    names = ['average', 'mid']
    
    wavelength = wavelength
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    date_title = '%s-%s-%s' % (year,month,day)
    
    vis = np.load('%s' % dataset)
    
    for i in range(2):
        
        fig = plt.figure(figsize=(15,9))
        ax = plt.gca()
        plt.title('SDO AIA %i.0 Angstrom %s [%s]' % (wavelength, date_title, titles[i]), y = 1.01, fontsize=25)
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
        plt.savefig('%s/%s_%i_visual_%s.jpeg' % (path_name, date, wavelength, names[i]))
        

def other(who, what, when):
    """This function does something.
 
    :param name: The name to use.
    :type name: str.
    :param state: Current state to be in.
    :type state: bool.
    :returns:  int -- the return code.
    :raises: AttributeError, KeyError
 
    """

    titles = ['Average', 'Middle-File']
    names = ['average', 'mid']
    
    wavelength = wavelength
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    date_title = '%s-%s-%s' % (year,month,day)
