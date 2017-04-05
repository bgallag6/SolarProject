#!/bin/bash

echo "This program loads in .FITS files, extracts the specified subregion, and returns the derotated datacube, exposure and time arrays." 

read -p "Enter a directory [ex. /media/solar/Gallagher]: " directory
read -p "Enter a date [ex. 20130626]: " date
read -p "Enter a wavelength [ex. 193]: " wavelength
read -p "Enter the xmin dimension [ex. -500]: " xmin
read -p "Enter the xmax dimension [ex. 500]: " xmax
read -p "Enter the ymin dimension [ex. -500]: " ymin
read -p "Enter the ymax dimension [ex. 500]: " ymax
read -p "Enter the number of processors [ex. 16]: " num

mpiexec -n $num python mpi_derotate.py $directory $date $wavelength $xmin $xmax $ymin $ymax

read -p "Enter the number of processors [ex. 16]: " num2