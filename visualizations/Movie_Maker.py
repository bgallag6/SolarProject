# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 20:56:05 2016

@author: Brendan
"""

import matplotlib.pyplot as plt

from sunpy.map import Map
import numpy as np
import scipy.signal
import matplotlib
#matplotlib.use('TkAgg') 	# NOTE: This is a MAC/OSX thing. Probably REMOVE for linux/Win
from matplotlib.widgets import  RectangleSelector
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from pylab import *
import glob
import sunpy
from sunpy.map import Map
from sunpy.image.coalignment import mapcube_coalign_by_match_template
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
from scipy.interpolate import interp1d
from scipy import signal
import numpy as np
import scipy.misc
import astropy.units as u
import h5py

#import matplotlib.animation as animation
import matplotlib.animation as manimation

# create a list of all the files. This is USER-DEFINED
flist = glob.glob('F:/Users/Brendan/Desktop/SolarProject/FITS/Movies/20130626/continuum/*.fits')
#flist = glob.glob('C:/Users/Brendan/Desktop/SDO/304/aia*.fits')
#flist = glob.glob('D:/1600part2/aia*.fits')
#flist2 = flist[0:25]
#nf = len(flist2)
nf = len(flist)


# Select the image that is the "middle" of our selection.
# We do this because the solar derotation algorithm operates centered on the 
# "middle" image
mid_file = np.int(np.floor(nf / 2))
m1 = Map(flist[mid_file])

# Get defined rectangle coords as integers
#x1=1400
#x2=1578
#y1=2400
#y2=2772

#"""
x1=2000
x2=3000
y1=2000
y2=3000
#"""

# Create an empty list
mc_list = []

# Use coordinates we just defined, extract the submaps from each AIA image, and store
# then in the empty list. This takes many minutes to complete.
print " "
print "Reading files and extracting submaps. This takes a while..."
print " "
for filename in flist:
   mc_list.append(Map(filename).submap([y1,y2]*u.pixel, [x1,x2]*u.pixel))  # ** SEE NOTE BELOW ***
   #mc_list.append(Map(filename).submap([y2,y1]*u.pixel, [x2,x1]*u.pixel))  # ** SEE NOTE BELOW ***

new_mapcube = Map(mc_list, cube=True)

dr = mapcube_solar_derotate(new_mapcube)
#dr = new_mapcube

"""
# full solar disk
temp = []

for i in range(1,64):
    temp[i+1]= flist[24*i]
# Select the image that is the "middle" of our selection.
# We do this because the solar derotation algorithm operates centered on the 
# "middle" image
mid_file = np.int(np.floor(nf / 2))
m1 = Map(flist[mid_file])


# Create an empty list
mc_list = []

# Use coordinates we just defined, extract the submaps from each AIA image, and store
# then in the empty list. This takes many minutes to complete.
print " "
print "Reading files and extracting submaps. This takes a while..."
print " "
for filename in temp:
   mc_list.append(Map(filename))  # ** SEE NOTE BELOW ***

new_mapcube = Map(mc_list, cube=True)

dr = mapcube_solar_derotate(new_mapcube)
"""


#cube = Map(mc_list, cube=True)   # doctest: +SKIP
#ani = cube.plot()   # doctest: +SKIP
#plt.show()   # doctest: +SKIP

#Plot the map at 1/2 original resolution

#cube = Map(files, cube=True)   # doctest: +SKIP
#ani = cube.plot(resample=[0.5, 0.5], colorbar=True)   # doctest: +SKIP
#plt.show()   # doctest: +SKIP

#Save an animation of the MapCube

#cube = Map(mc_list, cube=True)   # doctest: +SKIP

"""
## use these
ani = dr.plot()   # doctest: +SKIP

Writer = animation.writers['ffmpeg']   # doctest: +SKIP
writer = Writer(fps=60, metadata=dict(artist='SunPy'), bitrate=5000000)   # doctest: +SKIP

ani.save('C:/Users/Brendan/Desktop/try1.mp4', writer=writer, bbox_inches='tight')   # doctest: +SKIP
"""


#####
#possible other
#####

my_dpi=192 # dpi for retina display (detect monitor dpi automatically here: http://www.infobyip.com/detectmonitordpi.php)
#fig = plt.figure(figsize=(xsize/my_dpi, ysize/my_dpi), dpi=my_dpi)
fig = plt.figure(figsize=(12,10), dpi=175)
#print "The figure is the following size (in pixels):",fig
ax = fig.add_subplot(1,1,1)
#ax.get_xaxis().set_ticks([])
#ax.get_yaxis().set_ticks([])
savefigdict = {'bbox_inches' : 'tight'}

#Writer = manimation.writers['ffmpeg']   # doctest: +SKIP
#writer = Writer(fps=60, metadata=dict(artist='SunPy'), bitrate=5000000)   # doctest: +SKIP

ims = []
for i in range(nf):
    #im = plt.imshow(dr[:,:,i],cmap='sdoaia171',origin='lower',vmin=0,vmax=3000,extent=[0,xsize,0,ysize], interpolation=None, animated=True)
    #im = plt.imshow(dr[i].data,cmap='sdoaia171',origin='lower',vmin=0,vmax=3000, interpolation=None, animated=True)
    im = plt.imshow((dr[i].data)/1e3,cmap='gray',origin='lower',vmin=10,vmax=70, interpolation=None, animated=True)
    ims.append([im])

ani = manimation.ArtistAnimation(fig, ims, interval=500, blit=True)

#ani.save('C:/Users/Brendan/Desktop/test1.mp4', savefig_kwargs=savefigdict, writer=writer, dpi=my_dpi)