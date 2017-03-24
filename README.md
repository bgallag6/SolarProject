# SolarProject

## An Overview of how to generate heatmaps for a region. ##

1. Data is downloaded from the VSO for the desired region.
   - so far we have been working with either six or twelve-hour timespans, at 
     12 or 24-second cadence depending on what is available per particular wavelength.

2. The .fits files are first read into python.  The desired region is then specified and extracted from each file, put into a datacube, and the datacube is then de-rotated.  (The region's dimensions must take into account the trimming effect of the de-rotation.)

3. The image data is extracted from the cube and normalized by dividing through by the exposure time of each image.  The time-range of the dataset is also extracted, and both are saved in an HDF5 file.  

4. The intensity value for each pixel is extracted, and the timeseries is split into two-hour segments.  Each of these segments are processed through the Fast Fourier Transform, and then are averaged together to reduce the noise of the spectra.

5. To further reduce the noise and allow for a better fit, the region is then looped through by 3x3 pixel box.  The nine spectra in each are geometrically averaged and the result is assigned to the central pixel in the box.  

6. Each pixel's spectra is then fit to two models: the first made from a power-law-with-tail, and the second including an additional Gaussian component.

7. The six parameters from the combined model are extracted from the fits, as well as the chi-squared statistic.  In addition, the combined model fit and the averaged spectra are saved - to check in case of any errors.  

8. For each parameter, a heatmap is generated over the full region.


***
***

# A brief description on the organization of this repository + how the scripts should be used. #

### Main Folder: ###

- **download_fits.sh** : download .FITS files

- **derotate_slice_mpi2.py** : create datacube + derotate slices

- **merge_derotate_chunks.py** : merge the slices saved in previous file

- **heatmap_spectra_tool.py** : GUI for investigating individual pixel's spectra + fits

- **SolSpecShell.sh** : takes derotated cube + computes power spectra + memory maps result for fitting + fits spectra


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
