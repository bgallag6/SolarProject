# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 01:19:23 2017

@author: Brendan
"""

import numpy as np
from pylab import imshow, show
from timeit import default_timer as timer
from numba import autojit, vectorize, float64
from numbapro import cuda, vectorize
from numba import *
import numpy as np

import numpy as np
import scipy.signal
from scipy.interpolate import interp1d
from scipy import signal
import scipy.misc
import astropy.units as u
from scipy import fftpack  # doesnt work in module when called here???
from astropy.convolution import convolve, Box1DKernel
from numpy.random import randn
from mpi4py import MPI
import matplotlib.pyplot as plt
import math

from scipy import fftpack    



#"""
#@vectorize(["float32(float32,float32,float32,float32)","float64(float64,float64,float64,float64)"],target = 'cpu')
#def cpu_m2(f, a, n, c):
#    return a*f**(-n) + c

"""
#@vectorize(["float32(float32,float32,float32,float32)","float64(float64,float64,float64,float64)"],target = 'gpu')
@cuda.jit('float64(float64, float64, float64, float64)', device=True, inline=True)
def gpu_m2(f, a, n, c):
    return a*f**(-n) + c
 
@vectorize("float64(float64, float64, float64, float64)", target = 'gpu')
def gpu_m2_ufn(f,a,n,c):
    return gpu_m2(f,a,n,c)
"""


#@vectorize([float64(float64,float64,float64,float64,float64,float64,float64)], target='gpu')
#def gpu_m2_ufn(f,a,n,c,p,beta,sigma):
@vectorize([float64(float64,float64,float64,float64,float64,float64,float64)])
def gpu_m2_ufn(f,a,n,c,p,beta,sigma):
    #return a*f**(-n) + c + p*np.exp(-0.5(np.log(f)-(beta))**2/sigma**2)
    return a*f**(-n) + c + p*math.exp(-(math.log(f)-beta)**2/2*sigma**2)
    
def func_temp(f,a,n,c,p,beta,sigma):
    res = gpu_m2_ufn(f,a,n,c,p,beta,sigma)
    return res
    
    

f = np.linspace((1./5000.),(1./25.),1000).astype(np.float64)
a = 10**-8
n = 2.
c = 10**-3
p = 10**-3.5
beta = -5.2
sigma = 0.15
s = a*f**(-n) + c + p*np.exp(-(np.log(f)-beta)**2/2*sigma**2)
ds = 0.1*s



#threads_per_block = 128

#r = gpu_m2_ufn(f,a,n,c,p,beta,sigma)


M1_low = [-0.002, 0.3, -0.01, -0.01, -6.4, 0.05]
M1_high = [0.002, 6., 0.01, 0.01, 4.6, 0.8]

t0 = timer()

#nlfit_l, nlpcov_l = scipy.optimize.curve_fit(gpu_m2_ufn, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')
nlfit_l, nlpcov_l = scipy.optimize.curve_fit(func_temp, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')
    
#r = gpu_m2_ufn(f,a,n,c,p,beta,sigma)
t2 = timer()
print t2-t0

#f = np.linspace((1./5000.),(1./25.),10000).astype(np.float64)
#s = 10**-8*f**(-2.)+10**3.5

#gpu_m2(f,10**-8,-2.,10**-4)
#"""

#"""   
#@jit(nopython=True)
#def gpu_m2(f,a,n,c):
#    z = a*f**(-n) + c
#    return x
""" 
@jit(nopython=True)
def GaussPowerBase(f2, A2, n2, C2):
    return A2*f2**-n2 + C2


f = np.linspace((1./5000.),(1./25.),10000).astype(np.float64)
s = 10**-8*f**(-2.)+10**3.5

# assign equal weights to all parts of the curve
df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
df2 = np.zeros_like(f)
df2[0:len(df)] = df
df2[len(df2)-1] = df2[len(df2)-2]
ds = df2

threadsperblock = 256
blockspergrid = (f.size + (threadsperblock-1)) // threadsperblock

x1 = np.random.random(100).astype(np.float64)
x2 = np.random.random(100).astype(np.float64)
x3 = np.random.random(100).astype(np.float64)
M1_low = [-0.002, 0.3, -0.01]
M1_high = [0.002, 6., 0.01]

t0 = timer()
for i in range(len(x1)):
    nlfit_l, nlpcov_l = scipy.optimize.curve_fit(GaussPowerBase, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')
t1 = timer()
print t1-t0
"""    
"""
threadsperblock = 256

f = np.random.random(100).astype(np.float64)
s = 10**-8*f**(-2.)+10**3.5
# assign equal weights to all parts of the curve
df = np.log10(f[1:len(f)]) - np.log10(f[0:len(f)-1])
df2 = np.zeros_like(f)
df2[0:len(df)] = df
df2[len(df2)-1] = df2[len(df2)-2]
ds = df2

x1 = np.random.random(100).astype(np.float64)
x2 = np.random.random(100).astype(np.float64)
x3 = np.random.random(100).astype(np.float64)
M1_low = [-0.002, 0.3, -0.01]
M1_high = [0.002, 6., 0.01]

#t0 = timer()
#res1 = cpu_m2(x0,x1,x2,x3)
#t1 = timer()

t2 = timer()
#res2 = gpu_m2(x0,x1,x2,x3)
nlfit_l, nlpcov_l = scipy.optimize.curve_fit(gpu_m2, f, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')
t3 = timer()


print t3-t2
"""

 
"""   
@jit(nopython=True)
def nan_compact(x):
    out = np.empty_like(x)
    out_index = 0
    for element in x:
        if not np.isnan(element):
            out[out_index] = element
            out_index += 1
        return out[:out_index]
        
a = np.random.uniform(size=10000)
a[a < 0.2] = np.nan
np.testing.assert_equal(nan_compact(a), a[~np.isnan(a)])
"""



"""
##########################
# non-accelerated version
##########################
"""
"""
def mandel(x, y, max_iters):

  c = complex(x, y)
  z = 0.0j
  for i in range(max_iters):
    z = z*z + c
    if (z.real*z.real + z.imag*z.imag) >= 4:
      return i

  return max_iters
  
def create_fractal(min_x, max_x, min_y, max_y, image, iters):
  height = image.shape[0]
  width = image.shape[1]

  pixel_size_x = (max_x - min_x) / width
  pixel_size_y = (max_y - min_y) / height
    
  for x in range(width):
    real = min_x + x * pixel_size_x
    for y in range(height):
      imag = min_y + y * pixel_size_y
      color = mandel(real, imag, iters)
      image[y, x] = color

image = np.zeros((1024, 1536), dtype = np.uint8)
start = timer()
create_fractal(-2.0, 1.0, -1.0, 1.0, image, 20) 
dt = timer() - start

print "Mandelbrot created in %f s" % dt
imshow(image)
show()
"""



"""
##########################
# accelerated-CPU version
##########################
"""
"""
@autojit
def mandel(x, y, max_iters):

  c = complex(x, y)
  z = 0.0j
  for i in range(max_iters):
    z = z*z + c
    if (z.real*z.real + z.imag*z.imag) >= 4:
      return i

  return max_iters

@autojit
def create_fractal(min_x, max_x, min_y, max_y, image, iters):
  height = image.shape[0]
  width = image.shape[1]

  pixel_size_x = (max_x - min_x) / width
  pixel_size_y = (max_y - min_y) / height
    
  for x in range(width):
    real = min_x + x * pixel_size_x
    for y in range(height):
      imag = min_y + y * pixel_size_y
      color = mandel(real, imag, iters)
      image[y, x] = color

image = np.zeros((1024, 1536), dtype = np.uint8)
start = timer()
create_fractal(-2.0, 1.0, -1.0, 1.0, image, 20) 
dt = timer() - start

print "Mandelbrot created in %f s" % dt
imshow(image)
show()
"""



"""
##########################
# accelerated-GPU version
##########################
"""
"""
@autojit
def mandel(x, y, max_iters):

  c = complex(x, y)
  z = 0.0j
  for i in range(max_iters):
    z = z*z + c
    if (z.real*z.real + z.imag*z.imag) >= 4:
      return i

  return max_iters

mandel_gpu = cuda.jit(restype=uint32, argtypes=[f8, f8, uint32], device=True)(mandel)


@cuda.jit(argtypes=[f8, f8, f8, f8, uint8[:,:], uint32])
def mandel_kernel(min_x, max_x, min_y, max_y, image, iters):
  height = image.shape[0]
  width = image.shape[1]

  pixel_size_x = (max_x - min_x) / width
  pixel_size_y = (max_y - min_y) / height

  startX, startY = cuda.grid(2)
  gridX = cuda.gridDim.x * cuda.blockDim.x;
  gridY = cuda.gridDim.y * cuda.blockDim.y;

  for x in range(startX, width, gridX):
    real = min_x + x * pixel_size_x
    for y in range(startY, height, gridY):
      imag = min_y + y * pixel_size_y 
      image[y, x] = mandel_gpu(real, imag, iters)

gimage = np.zeros((1024, 1536), dtype = np.uint8)
blockdim = (32, 8)
griddim = (32,16)

start = timer()
d_image = cuda.to_device(gimage)
mandel_kernel[griddim, blockdim](-2.0, 1.0, -1.0, 1.0, d_image, 20) 
d_image.to_host()
dt = timer() - start

print "Mandelbrot created on GPU in %f s" % dt

imshow(gimage)
show()
"""