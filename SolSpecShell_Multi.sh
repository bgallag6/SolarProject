#!/bin/bash

echo "The process of this program:
1) Load in .FITS files from desired wavelength
2) Creates datacube and applies derotation algorithm
3) Extracts pixel intensity values, and exposure duration from each .FITS image
4) Power-spectra are computed from extracted timeseries via use of the Fast Fourier Transform
5) A memory-mapped copy of the power-spectra array is created in order to be passed to MPI
6) Two models, M1 and M2, are fitted to the spectra and the model parameters are extracted
7) Heatmaps and histograms are generated from extracted model parameters" 

read -p "Enter a directory [ex. /media/solar/Gallagher]: " directory
read -p "Enter a date [ex. 20130626]: " date
read -p "Enter how many wavelengths [ex. 5]: " num_wave

if [ $num_wave -eq 1 ]
then 
	read -p "Wavelength 1?: " wave1
	echo "wavelengths: $wave1"

elif [ $num_wave -eq 2 ]
then 
	read -p "Wavelength 1?: " wave1
	read -p "Wavelength 2?: " wave2
	echo "wavelengths: $wave1 $wave2"

elif [ $num_wave -eq 3 ]
then 
	read -p "Wavelength 1?: " wave1
	read -p "Wavelength 2?: " wave2
	read -p "Wavelength 3?: " wave3
	echo "wavelengths: $wave1 $wave2 $wave3"

elif [ $num_wave -eq 4 ]
then 
	read -p "Wavelength 1?: " wave1
	read -p "Wavelength 2?: " wave2
	read -p "Wavelength 3?: " wave3
	read -p "Wavelength 4?: " wave4
	echo "wavelengths: $wave1 $wave2 $wave3 $wave4"

elif [ $num_wave -eq 5 ]
then 
	read -p "Wavelength 1?: " wave1
	read -p "Wavelength 2?: " wave2
	read -p "Wavelength 3?: " wave3
	read -p "Wavelength 4?: " wave4
	read -p "Wavelength 5?: " wave5
	echo "wavelengths: $wave1 $wave2 $wave3 $wave4 $wave5"
fi

read -p "Enter the number of processors [ex. 16]: " num

if [ $num_wave -eq 1 ]
then 
	python SolSpec_Call.py $directory $date $wave1

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave1

elif [ $num_wave -eq 2 ]
then 
	python SolSpec_Call.py $directory $date $wave1

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave1

	python SolSpec_Call.py $directory $date $wave2

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave2

elif [ $num_wave -eq 3 ]
then 
	python SolSpec_Call.py $directory $date $wave1

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave1

	python SolSpec_Call.py $directory $date $wave2

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave2

	python SolSpec_Call.py $directory $date $wave3

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave3

elif [ $num_wave -eq 4 ]
then 
	python SolSpec_Call.py $directory $date $wave1

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave1

	python SolSpec_Call.py $directory $date $wave2

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave2

	python SolSpec_Call.py $directory $date $wave3

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave3

	python SolSpec_Call.py $directory $date $wave4

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave4

elif [ $num_wave -eq 5 ]
then 
	python SolSpec_Call.py $directory $date $wave1

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave1

	python SolSpec_Call.py $directory $date $wave2

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave2

	python SolSpec_Call.py $directory $date $wave3

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave3

	python SolSpec_Call.py $directory $date $wave4

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave4

	python SolSpec_Call.py $directory $date $wave5

	mpiexec -n $num python Spec_fit_mpi.py $directory $date $wave5
fi