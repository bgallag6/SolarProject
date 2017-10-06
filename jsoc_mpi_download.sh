#!/bin/bash

echo "The process of this program:
1) Scrubs file download links from JSOC url 
2) Splits download requests among processors
3) Renames files to match Python/VSO downloads
4) Saves downloaded files to correct folder

Folder structure is currently: 
*directory*/FITS/*date*/*wavelength*/{}.fits
"

read -p "Enter the JSOC url: " url_jsoc
read -p "Enter a directory for files to be saved [ex. /media/solar/Gallagher]: " directory
read -p "Enter the date [ex. 20130626]: " date
read -p "Enter the number of processors [ex. 16]: " num

:: mpiexec -n $num python mpi_jsoc_download.py $url_jsoc $directory $date

:: mpiexec -n $num python mpi_jsoc_euv_download.py $url_jsoc $directory $date

:: mpiexec -n $num python mpi_jsoc_hmi_download.py $url_jsoc $directory $date

:: mpiexec -n $num python mpi_jsoc_continuum_download.py $url_jsoc $directory $date

mpiexec -n $num python mpi_jsoc_download_complete.py $url_jsoc $directory $date