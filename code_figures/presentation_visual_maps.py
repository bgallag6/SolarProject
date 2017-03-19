# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 08:56:36 2017

@author: Brendan
"""

import matplotlib.pyplot as plt
from matplotlib import patches
import astropy.units as u
import sunpy.map



# Define a region of interest
length = 479 * u.arcsec
x0 =  0 * u.arcsec
y0 = 0 * u.arcsec


# Create a SunPy Map, and a second submap over the region of interest.
smap = sunpy.map.Map('F:/Users/Brendan/Desktop/SolarProject/FITS/20130530/193/aia_lev1_193a_2013_05_30t00_00_06_84z_image_lev1.fits')
submap = smap.submap(u.Quantity([x0 - length, x0 + length]),
                     u.Quantity([y0 - length, y0 + length]))


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