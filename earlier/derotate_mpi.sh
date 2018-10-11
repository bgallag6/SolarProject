#!/bin/bash

echo "Program Steps: 
1) Use MPI to differentially derotate FITS images
2) Save chunks of derotated FITS images as numpy files
3) Load the separate chunk files, merge them and save the final datacube
4) Also save the exposure and time arrays, as well as normalized averaged visual image
"""

read -p "Enter a directory [ex. /media/solar/Gallagher]: " directory
read -p "Enter a date [ex. 20130626]: " date
read -p "Enter a wavelength [ex. 193]: " wavelength
read -p "Enter x1 coordinate [ex. -500]: " x1
read -p "Enter x2 coordinate [ex. 500]: " x2
read -p "Enter y1 coordinate [ex. -500]: " y1
read -p "Enter y2 coordinate [ex. 500]: " y2
read -p "Enter the number of processors [ex. 16]: " num

mpiexec -n $num python diff_derotate_mpi.py $directory $date $wavelength $x1 $x2 $y1 $y2

python merge_diff_chunks.py $directory $date $wavelength $num