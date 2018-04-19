#!/bin/bash

echo "The process of this program:
1) Load in .FITS files from desired wavelength
2) Creates datacube and applies derotation algorithm
3) Extracts pixel intensity values, and exposure duration from each .FITS image
4) A memory-mapped copy of the derotated array is created in order to be passed to MPI
5) Power-spectra are computed from extracted timeseries via use of the Fast Fourier Transform
6) A memory-mapped copy of the power-spectra array is created in order to be passed to MPI
7) Two models, M1 and M2, are fitted to the spectra and the model parameters are extracted"

read -p "Enter the number of processors [ex. 16]: " num

:: mpiexec -n $num python part1B_diff_derotate.py

:: python part15B_merge_diff_chunks.py $num

:: mpiexec -n $num python part2B_fft_avg_mpi_Gather.py

mpiexec -n $num python part3B_spec_fit_mpi_master.py
