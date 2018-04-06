# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 08:56:36 2017

@author: Brendan
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import astropy.units as u
import sunpy.map
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

plt.rcParams["font.family"] = "Times New Roman"

# Define a region of interest
length = 479 * u.arcsec
x0 =  0 * u.arcsec
y0 = 0 * u.arcsec

sq = 1000 * u.arcsec


# Create a SunPy Map, and a second submap over the region of interest.
smap = sunpy.map.Map('F:/FITS/20130626/171/aia_lev1_171a_2013_06_26t00_00_11_34z_image_lev1.fits').submap(u.Quantity([-sq,sq]),u.Quantity([-sq,sq]))
submap = smap.submap(u.Quantity([x0 - length, x0 + length]),
                     u.Quantity([y0 - length, y0 + length]))
            
data = smap.data
ex = (smap.exposure_time).value


"""
# Create a new matplotlib figure, larger than default.
fig = plt.figure(figsize=(22,10))

# Add a first Axis, using the WCS from the map.
ax1 = plt.subplot2grid((1,11),(0, 0), colspan=5, rowspan=1, projection=smap)
#ax1 = fig.add_subplot(2,1,1, projection=smap)

# Plot the Map on the axes with default settings.
smap.plot()
ax1.set_title(r'SDO AIA 193$\AA$ 2013/05/30 00:00:06 UTC', fontsize=20, y=1.01)

# Define a region to highlight with a box
# We have to convert the region of interest to degress, and then get the raw values.
bottom_left = u.Quantity([x0 - length, y0 - length])
length2 = length * 2

# Draw a box on the image
smap.draw_rectangle(bottom_left, length2, length2)

# Create a second axis on the plot.
ax2 = plt.subplot2grid((1,11),(0, 6), colspan=5, rowspan=1, projection=submap)
#ax2 = fig.add_subplot(2,2,1, projection=submap)

submap.plot()

# Add a overlay grid.
submap.draw_grid(grid_spacing=10*u.deg)

# Change the title.
ax2.set_title('Zoomed View - 1600x1600-Pixel Region', fontsize=20, y=1.01)


plt.show()
#plt.savefig('C:/Users/Brendan/Desktop/193_maps_full.pdf', format='pdf')
#plt.savefig('C:/Users/Brendan/Desktop/193_maps_full.jpeg')
"""

font_size = 27

#x_ticks = [-1000,-500,0,500,1000]
#y_ticks = [-1000,-500,0,500,1000]
x_ticks = [-958,-479,0,479,958]
y_ticks = [-958,-479,0,479,958]
x_ind = [-1600,-800,0,800,1600]
y_ind = [-1600,-800,0,800,1600]
c_ticks = [100,200,300,400,500,600,700,800]

trim = (4096-1600)/2

#data_ex_norm = data/ex
data_ex_norm = data
data_trim = data_ex_norm[trim:-trim,trim:-trim]
v_min = np.percentile(data_trim,1.)
v_max = np.percentile(data_trim,99.)

fig = plt.figure(figsize=(12,10))
# Add a first Axis, using the WCS from the map.
ax1 = plt.gca()
#ax1 = plt.subplot2grid((1,1),(0, 0), colspan=1, rowspan=1, projection=smap)
#ax1 = fig.add_subplot(2,1,1, projection=smap)
# Plot the Map on the axes with default settings.
#norm = colors.Normalize(vmin=smap.min(), vmax=smap.mean() + 3 *smap.std())
norm = colors.Normalize(vmin=v_min, vmax=v_max)
im = smap.plot(norm=norm)

#im = smap.plot(cmap = 'sdoaia171', vmin = 100, vmax = 1800)
#ax1.set_title(r'SDO AIA 171$\AA$ 2013/06/26 05:59:59 UTC', fontsize=font_size, y=1.01)
ax1.set_title(r'171 $\AA$ - 05:59:59 UT', fontsize=font_size, y=1.02, fontname="Times New Roman")
ax1.set_xlabel('', visible=False, fontname="Times New Roman")
ax1.set_ylabel('', visible=False, fontname="Times New Roman")
#plt.xticks(x_ticks,fontsize=font_size)
#plt.yticks(y_ticks,fontsize=font_size)
plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
ax1.tick_params(axis='both', which='major', pad=10)
#ax1.tick_params(axis='both', which='both', labelleft='off', labelbottom='off')
# Define a region to highlight with a box
# We have to convert the region of interest to degress, and then get the raw values.
bottom_left = u.Quantity([x0 - length, y0 - length])
length2 = length * 2
# Draw a box on the image
smap.draw_rectangle(bottom_left, length2, length2, linewidth=3.)
#plt.colorbar()
#divider = make_axes_locatable(ax1)
#cax = divider.append_axes("right", size="3%", pad=0.07)
#cbar = plt.colorbar(im,cax=cax)
#cbar.set_label('Intensity', size=20, labelpad=10)
#cbar.ax.tick_params(labelsize=font_size, pad=5)
#plt.savefig('C:/Users/Brendan/Desktop/171_full_disc_arcsecG.pdf', format='pdf')