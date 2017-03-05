#!/bin/bash

echo "The process of this program:
1) Load in .FITS files from desired wavelength
2) Creates datacube and applies derotation algorithm
3) Extracts pixel intensity values, and exposure duration from each .FITS image
4) Power-spectra are computed from extracted timeseries via use of the Fast Fourier Transform
5) A memory-mapped copy of the power-spectra array is created in order to be passed to MPI
6) Two models, M1 and M2, are fitted to the spectra and the model parameters are extracted
7) Heatmaps and histograms are generated from extracted model parameters 

read -p "Enter a directory: " directory
read -p "Enter a date: " date
read -p "Enter a wavelength: " wavelength
read -p "Enter the number of processors: " num

python SolSpec_Call.py $directory $date $wavelength

mpiexec -n $num python Spec_fit_mpi.py $directory $date $wavelength

python SolSpec_Call_Heatmaps.py $directory $date $wavelength
