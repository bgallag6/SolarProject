#!/bin/bash

echo "This program searches the provided directory for the .FITS files already 
downloaded, creates an array containing the images still needed, 
then downloaded those to the same directory."

read -p "Enter a directory [ex. /media/solar/Gallagher]: " directory
read -p "Enter a date without slashes [ex. 20130530]: " date
read -p "Enter a start time [ex. 00:00:00]: " start_time
read -p "Enter an end time [ex. 11:59:59]: " end_time
read -p "Enter a wavelength [ex. 193]: " wavelength
read -p "Enter the number of processors [ex. 16]: " num

python create_arr_need.py $directory $date $start_time $end_time $wavelength

mpiexec -n $num python mpi2_download_fits.py $directory $date $wavelength