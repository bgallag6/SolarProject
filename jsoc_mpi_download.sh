#!/bin/bash

echo "The process of this program:
1) Scrubs file download links from JSOC url 
2) Splits download requests among processors
3) Renames files to match Python/VSO downloads
4) Saves downloaded files to correct folder

Folder structure is currently: 
*directory*/FITS/*date*/*wavelength*/{}.fits
"

read -p "Enter a directory for files to be saved [ex. /media/solar/Gallagher]: " directory
read -p "Enter the date [ex. 20130626]: " date
read -p "Enter the wavelength [ex. 1600]: " wavelength
read -p "Enter the start time [ex. 00:00]: " tstart
read -p "Enter the duration in hours [ex. 12]: " duration
read -p "Enter the number of processors [ex. 16]: " num

python make_fits_directory.py $directory $date $wavelength

python -u JSOC_request_url.py $directory $date $wavelength $tstart $duration

mpiexec -n $num python -u mpi_jsoc_download_complete.py $directory $date $wavelength