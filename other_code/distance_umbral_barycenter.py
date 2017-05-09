# -*- coding: utf-8 -*-
"""
Created on Mon May 01 13:27:14 2017

@author: Brendan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import astropy.units as u
import sunpy
from sunpy.map import Map
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm

#H = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20130626/PCB/1600/param.npy')
H = np.load('F:/Users/Brendan/Desktop/SolarProject/DATA/Output/20111210/1600/param.npy')

#L = H[4,180:310,250:400]  # 20130626 PCB 1600 (76x, 59y)
L = H[4,92:163,90:167]  # 20111210 (38x, 37y)
L = (1./np.exp(L))/60.

"""
y = np.zeros((7,7))
z = np.zeros((7,7))

for i in range(7):
    for j in range(7):
        posx = np.abs(3-j)
        posy = np.abs(3-i)
        pos = np.max([posx,posy])
        if pos == 0:
            y[i][j] = 10
        elif pos == 1:
            y[i][j] = 5
        elif pos == 2:
            y[i][j] = 1
        if (posx**2 + posy**2) < 1**2:
            z[i][j] = 10
        elif (posx**2 + posy**2) < 2**2 and (posx**2 + posy**2) >= 1**2:
            z[i][j] = 5
        elif (posx**2 + posy**2) >= 2**2 and (posx**2 + posy**2) < 3**2:
            z[i][j] = 2
        elif (posx**2 + posy**2) >= 3**2 and (posx**2 + posy**2) < 4**2:
            z[i][j] = 1
            
#plt.figure()            
#plt.imshow(y)

#plt.figure()
#plt.imshow(z)    
#z = ((5-2)**2 + (5-2)**2)**(.5)

z1 = 0
z2 = 0
z3 = 0
z1c = 0
z2c = 0
z3c = 0

for i in range(7):
    for j in range(7):
        posx = np.abs(3-j)
        posy = np.abs(3-i)
        pos = np.max([posx,posy])
        if (posx**2 + posy**2) < 2**2 and (posx**2 + posy**2) >= 1**2:
            z1 += z[i][j]
            z1c += 1
        elif (posx**2 + posy**2) >= 2**2 and (posx**2 + posy**2) < 3**2:
            z2 += z[i][j]
            z2c += 1
        elif (posx**2 + posy**2) >= 3**2 and (posx**2 + posy**2) < 4**2:
            z3 += z[i][j]
            z3c += 1
            
z1 = z1/z1c
z2 = z2/z2c
z3 = z3/z3c

#x1 = np.linspace(0,L.shape[1],L.shape[1])
#x2 = np.linspace(0,L.shape[0],L.shape[0])
#x1,x2 = np.meshgrid(x1,x2)
#R = (x1-76)**2 + (x2-76)**2
"""

pix = 25

#x_cen = 76 # 20130626
#y_cen = 59

x_cen = 38 # 20111210
y_cen = 37

z = np.zeros((pix))
zcount = np.zeros((pix))


for i in range(L.shape[0]):
    for j in range(L.shape[1]):
        
        posx = np.abs(x_cen-j)
        posy = np.abs(y_cen-i)
        pos = np.max([posx,posy])
        
        for c in range(len(z)):
            if (posx**2 + posy**2) < (c+1)**2 and (posx**2 + posy**2) >= c**2:
                z[c] += L[i][j]
                zcount[c] += 1


z = z/zcount

print z

lw = 2.

#r1 = 10 # 20130626
#r2 = 20

r1 = 7 # 20111210
r2 = 15

plt.figure(figsize=(20,7.2))
ax1 = plt.subplot2grid((1,11),(0, 0), colspan=5, rowspan=1)
#ax1.set_title(r'20130626 1600$\AA$: Period [min]', y=1.01, fontsize=17)
ax1.set_title(r'20111210 1600$\AA$: Period [min]', y=1.01, fontsize=17)
ax = plt.gca()  # get current axis -- to set colorbar 
#plt.title(r'%s: %i $\AA$  [%s]' % (date_title, wavelength, titles[i]), y = 1.01, fontsize=25)
im = ax.imshow(L, cmap='jet_r')
circle = plt.Circle((x_cen, y_cen), r1, linewidth=lw, color='white', fill=False, linestyle='dashed')
circleB = plt.Circle((x_cen, y_cen), r1, linewidth=lw+2.5, color='black', fill=False, linestyle='dashed')
circle2 = plt.Circle((x_cen, y_cen), r2, linewidth=lw, color='white', fill=False, linestyle='dashed')
circle2B = plt.Circle((x_cen, y_cen), r2, linewidth=lw+2.5, color='black', fill=False, linestyle='dashed')
ax.add_artist(circleB)
ax.add_artist(circle)
ax.add_artist(circle2B)
ax.add_artist(circle2)
divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
cax = divider.append_axes("right", size="3%", pad=0.07)
if i == 0:
    cbar = plt.colorbar(im,cax=cax, format='%0.2e')
else:
    cbar = plt.colorbar(im,cax=cax)
#cbar.set_label('%s' % cbar_labels[i], size=20, labelpad=10)
cbar.ax.tick_params(labelsize=17, pad=5) 
        
ax2 = plt.subplot2grid((1,11),(0, 6), colspan=5, rowspan=1)
ax2.set_title('Period [min] vs Distance from Umbral Barycenter', y=1.01)
plt.vlines(r1, 0,z.max(), linestyle='dashed', label='r1')
plt.vlines(r2, 0,z.max(), linestyle='dotted', label='r2')
plt.xlabel('Pixels', fontsize=13)
plt.ylabel('Period [min]', fontsize=13)
plt.legend(loc='upper left')
plt.plot(z)

#plt.savefig('C:/Users/Brendan/Desktop/20111210_1600_distance_period.pdf', format='pdf')
