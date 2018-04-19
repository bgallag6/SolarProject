# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 11:21:00 2018

@author: Brendan
"""


import numpy as np
from mpi4py import MPI
                 

def squareN( x ):
    return x+1
    
DIETAG = 9999 

class Work():
    def __init__(self, work_items):
        self.work_items = work_items[:] 
 
    def get_next_item(self):
        if len(self.work_items) == 0:
            return None
        return self.work_items.pop()
 
def master(wi):
    all_data = np.zeros((len(wi)))
    size = MPI.COMM_WORLD.Get_size()
    current_work = Work(wi) 
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    for i in range(1, size): 
        anext = current_work.get_next_item() 
        if not anext: break
        comm.send(obj=anext, dest=i, tag=anext)
 
    while 1:
        anext = current_work.get_next_item()
        print("next = ", anext, flush=True)
        if anext == None: break
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        #all_data.append(data)
        tg = status.Get_tag()
        print("node: 0  received tag: ", tg, flush=True)
        all_data[tg] = data
        comm.send(obj=anext, dest=status.Get_source(), tag=anext)
        print("node: 0  sent tag: ", anext, flush=True)
 
    for i in range(1,size):
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        #all_data.append(data)
        tg = status.Get_tag()
        print("node: 0  received tag: ", tg, flush=True)
        all_data[tg] = data
    
    for i in range(1,size):
        comm.send(obj=None, dest=i, tag=DIETAG)
     
    return all_data
        
    
def slave():
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    while 1:
        data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        rnk = comm.Get_rank()
        print("node: ", rnk, " tag: ", status.Get_tag(), flush=True)
        if status.Get_tag() == DIETAG: break
        y = squareN(data)
        comm.send(obj=y, dest=0, tag=status.Get_tag())
    

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size() 

if rank == 0:
    test_dat = [int(j) for j in range(10)]
    print(test_dat, flush=True)
    all_dat = master(test_dat)
    print(all_dat, flush=True)
else:
    slave()