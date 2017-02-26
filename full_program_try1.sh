#!/bin/bash

python SolSpec_Call.py

mpiexec -n 4 python Spec_fit_mpi.py

python SolSpec_Call_Heatmaps.py
