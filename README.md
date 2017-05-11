# SolarProject

## An Overview of how to generate heatmaps for a region. ##

1. Data is downloaded from the VSO for the desired region.
   - so far we have been working with either six or twelve-hour timespans, at 
     12 or 24-second cadence depending on what is available per particular wavelength.
   - the wavelengths we have investigated are: 171A, 193A, 211A, 304A, 1600A

2. The .FITS files are first read into python.  The desired region is then specified and extracted from each file, put into a datacube, and the datacube is then de-rotated.  (The region's dimensions must take into account the trimming effect of the de-rotation.)
   - the updated version of this step involves splitting the cube into slices and derotating each separately, 
     so that the correct derotation shifts are applied to each slice, and then merging the slices back into the full region

3. The image data is extracted from the cube and normalized by the exposure time of each image.  The time-range of the dataset is also extracted.  

4. The intensity value for each pixel is extracted, and the timeseries is split into two-hour segments.  For each of these segments the Fast Fourier Transform is computed.  The resulting power spectra are averaged together to reduce the noise level.

5. To further reduce the noise and allow for more efficient fitting, the region is then looped through by 3x3 pixel box.  The nine spectra in each are averaged and the result is assigned to the central pixel in the box.  

6. Each pixel's spectra is then fit to two models: the first a power law with tail, the second including an additional Gaussian component.

7. The six parameters from the combined model are extracted from the fits.  The f-statistic is calculated from the f-test of two models and their respective degrees of freedom.

8. For each parameter, a heatmap is generated for the full region.

9. A signicance test is applied to the Gaussian component parameters, masking the pixels whose p-value is above a designated threshold.  This step 'filters' the pixels whose spectra do not contain a 'significant' Gaussian component.  


***
***

# A brief description on the organization of this repository + how the scripts should be used. #

### Main Folder: ###

- **download_fits.sh** : download .FITS files

- **derotate_slice_mpi2.py** : create datacube + derotate slices

- **merge_derotate_chunks.py** : merge the slices saved in previous file

- **heatmap_spectra_tool.py** : GUI for investigating individual pixel's spectra + fits

- **SolSpecShell.sh** : takes derotated cube + computes power spectra + memory maps result for fitting + fits spectra

- **SolSpecShell_Multi.sh** : same as SolSpecShell.sh, but allows for specification of multiple wavelengths to process - passes them through all steps sequentially

- **create_arr_need.py** + **mpi2_download_fits.py** : uses MPI for multithreaded downloaded requests

- **FFT_avg_mpi.py** : MPI implementation of FFT-computing of power spectra (not currently in use)

- **SolSpec.py** : module containing functions: heatmap generation, derotation, FFT-computation, and memory-map creation of spectra cube

- **SolSpec_Call.py** : script calling the latter 3 functions in SolSpec.py

- **SolSpec_Call_Heatmaps.py** : script calling the heatmap generation function in SolSpec.py

- **Spec_fit_mpi.py** : spectra-fitting routine, called by SolSpecShell.sh

- **Spec_fit_mpi_4d.py** : spectra-fitting routine for time-evolution of regions, called by SolSpecShell.sh




### 'code_figures' Folder: ###

- scripts generating figures for the paper / presentation


### 'other_code' Folder: ###

- random other scripts, including older programs and ones I'm working on


### 'testing' Folder: ###

- scripts used to test various parts of the overall program


### 'validation' Folder: ###

- scripts generating spectra + fits to be cross-validated 


### 'visualizations' Folder: ###

- scripts generating some type of visualization used to analyze data
