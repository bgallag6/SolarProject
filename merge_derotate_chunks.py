# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 18:36:12 2017

@author: Brendan
"""

import numpy as np

#size = 16
size = 4

newData_p = []
for i in range(size):
    #temp = np.load('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/chunk_%i_of_16.npy' % i)
    temp = np.load('C:/Users/Brendan/Desktop/derotate_save/chunk_%i_of_4.npy' % i)
    newData_p.append(temp)


xmin = int(np.min([newData_p[r].shape[2] for r in range(size)]))  # this shouldn't be a huge difference - usually 2-3 pixels separate max/min
#for r in range(len(newData_p)):
    #temp_max[r] = newData_p[r].shape[2]
#xmin = int(np.min(temp_max))
#for s in range(len(newData_p)):
for s in range(size):
    newData_p[s] = newData_p[s][:,:,0:xmin]   # ugh, the trimming should probably be off of both sides of the array
num_f = newData_p[0].shape[0]
#y_len = newData_p[0].shape[1]
y_min = int(np.min([newData_p[w].shape[1] for w in range(size)]))  # this shouldn't be a huge difference - usually 2-3 pixels separate max/min
#y_len = 0
#for q in range(1,(size/2)):
#    y_len += np.min(y_min[(2*q)-1:(2*q)])
#y_full_dim
y_len = y_min
#y_len = int(np.min([newData_p[r].shape[1] for r in range(size)]))
#overlap = np.zeros((num_f,20,xmin))
overlap = int(np.floor(12*1.67))
#full_arr = np.zeros((num_f,(len(newData_p)*y_len)-(overlap*(len(newData_p)-1)),xmin))
full_arr = np.zeros((num_f,(size*y_len)-(overlap*(size-1)),xmin))

# maybe find coordinates that lead to even dimensions when dividing by size

weight_avg_top = np.array([np.linspace(20,1,20) for i in range(1666)])
weight_avg_bot = np.array([np.linspace(1,20,20) for i in range(1666)])
weight_avg_top = np.transpose(weight_avg_top)
weight_avg_bot = np.transpose(weight_avg_bot)
weight_avg_top = np.array([weight_avg_top for i in range(13)])
weight_avg_bot = np.array([weight_avg_bot for i in range(13)])

#for t in range(len(newData_p)-1):
for t in range(size-1):
    if t == 0:
        full_arr[:,0:(t+1)*(y_len-overlap),:] = newData_p[t][:,0:y_len-overlap,:]
        full_arr[:,(t+1)*(y_len-overlap):(t+1)*(y_len)-(t*overlap),:] = (newData_p[t][:,y_len-overlap:y_len,:]*weight_avg_top + newData_p[t+1][:,0:overlap,:]*weight_avg_bot) / 21.
    else:
        full_arr[:,(t*y_len)-((t-1)*overlap):(t+1)*(y_len-overlap),:] = newData_p[t][:,overlap:y_len-overlap,:]
        full_arr[:,(t+1)*(y_len-overlap):(t+1)*(y_len)-(t*overlap),:] = (newData_p[t][:,y_len-overlap:y_len,:]*weight_avg_top + newData_p[t+1][:,0:overlap,:]*weight_avg_bot) / 21.
full_arr[:,(t+1)*y_len-(t*overlap):(t+2)*y_len-((t+1)*overlap), :] = newData_p[size-1][:,overlap:y_len,:] 
    
#stack_p = np.hstack(newData_p)
#print stack_p.shape
#print "Derotated cube shape: (%s,%s,%s)" % full_arr.shape
#np.save('C:/Users/Brendan/Desktop/derotate_mpi_8seg_F.npy',full_arr)
#np.save('C:/Users/Brendan/Desktop/time.npy', newData_t[0])
#np.save('C:/Users/Brendan/Desktop/exposure.npy', newData_e[0])
#np.save('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/derotated.npy', full_arr)
np.save('C:/Users/Brendan/Desktop/derotate_save/derotated_weight.npy', full_arr)
#np.save('%s/DATA/Temp/%s/%i/time.npy' % (directory, date, wavelength), newData_t[0])
#np.save('%s/DATA/Temp/%s/%i/exposure.npy' % (directory, date, wavelength), newData_e[0])
#"""
# calculate the average-intensity image of the timeseries 
AVG = np.average(full_arr,axis=0)  #hmm, this is not normalized - doesn't really matter I don't think

# determine the middle file of the timeseries
mid_num = (full_arr.shape[0]/2)
mid = full_arr[mid_num]

print "Middle file is number %i" % mid_num

print "(100,100): %i ~ %i" % (AVG[100][100], mid[100][100])  # check values are reasonably close

# check the two image sizes agree
print " Average Image Dimensions = %i, %i" % (AVG.shape[0], AVG.shape[1])
print " Middle Image Dimensions = %i, %i" % (mid.shape[0], mid.shape[1])

# store average and middle images in array
visual = np.zeros((2,AVG.shape[0],AVG.shape[1]))
visual[0] = AVG
visual[1] = mid

print visual.shape  # check array size agrees with expected

# save visual-image array
#np.save('C:/Users/Brendan/Desktop/visual.npy', visual)
#np.save('/mnt/data-solar/Gallagher/DATA/Temp/20130626/171/visual.npy', visual)
np.save('C:/Users/Brendan/Desktop/derotate_save/visual_weight.npy', visual)
#np.save('%s/DATA/Output/%s/%i/visual.npy' % (directory, date, wavelength), visual)