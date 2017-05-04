## Solar Project Repository Changelog ##

4/6 : 
- added r-value calculation to Spec_fit_mpi.py and Spec_fit_mpi_4d.py
- increased parameter array size to 9 to hold r-value
- changed Spec_fit_mpi.py final time printout format from seconds to hours:minutes:seconds
- changed a bunch of "for i in range(0,...)" to "for i in range(...)"

4/7 : 
- added mean, median, +/- standard deviations statistics to histograms

4/8 :
- added M1 curve fitting + parameter display to fft_portion_select.py
- changed fft_portion_select from sample data to any data

- commented-out trimming timeseries if any negative intensity values (led to empty arrays - wasn't necessary?)

4/9 :
- added 3x3 pixel-box averaging spectra example + 3x3 pixel box visual example to Example_process_6_segments_FFT.py

- added 3x3 pixel-box averaging to fft_overlap in SolSpec.py

4/10 :
- changed heatmaps_4d.py colorscale to continuous vs discrete

4/12 : 
- added final printout of region and wavelength to Spec_fit_mpi.py and Spec_fit_mpi_4d.py
- commented-out uncertainty calculation in Spec_fit_mpi.py and Spec_fit_mpi_4d.py

- removed the prompt for the date with and without slashes from download_fits.sh
- replaced date-slash variable in create_arr_need.py and mpi2_download_fits.py with string creation with slashes from regular date

4/13 : 
- added r-value plot w/ stats to heatmap routine in SolSpec.py  (older parameter arrays won't work with this - since the r-values were extracted until recently)

4/17 : 
- indented the printout of the program time and saving of the parameter array in Spec_fit_mpi.py
(I guess it was trying to save the stacked array 16 times, and it wasn't defined for 15/16 processors.)

4/18 : 
- commented out upper/lower sigma + added 'mode' excluding 1st/last bins to histograms in SolSpec.py

4/28 : 
- changed the method of extracting time from FITS files in SolSpec.py

4/30:
- corrected time-estimation method in Spec_fit_mpi_4d.py
- added list sort function to glob list of FITS images in SolSpec.py

