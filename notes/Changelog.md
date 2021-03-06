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

5/5:
- added center weighted averaging to 3x3 pixel averaging in SolSpec.py

5/8:
- added p-value inversion mask contour overlay to visual images in SolSpec.py
- added sunspot umbra contour line to SolSpec.py

5/9:
- indented last 10 lines of Spec_fit_mpi.py - so only rank 0 is trying to save param.npy

5/15:
- for test: changed "for i in range(0,1)" - 1 2hour segment (still using 6 segments in call for freqs)
  also changed all the 3x3 averaging stuff back from 5x5... array allocating...

5/19:
- added calculating the standard deviation of the 3x3 pixel box's spectra as a function of frequency, to use as uncertainties 'ds' in spectra fitting [SolSpec.py]
- added option to constrain parameter bounds of 2nd optimization, since initial parameter guesses are specified [Spec_fit_mpi.py]
- changed a bunch of formatting of heatmaps [SolSpec.py]
- added creation of memory map of standard deviation array [SolSpec.py]
- added loading / splitting into chunks and passing to function of standard deviation memory mapped array, and using as each pixels 'ds' [Spec_fit_mpi.py]

8/10:
- removed a bunch of module imports + deleted unused / replaced lines [heatmap_spectra_tool.py]

8/12:
- added "bbox_inches='tight' to figure plots (removes whitespace around borders) - use scale vs. width in latex files [SolSpec.py]
- added rollover calculation to spec_fit.py 
- included r-value plot + histogram in loop, separated p-value plot creation from loop: heatmaps script
- included rollover plot + histogram in loop - for those parameter arrays that have it: heatmaps script
- moved umbra/PPV overlay into if-statement so no global error if not defined: heatmaps script

8/22:
- added Gaussian fit of Gaussian location to SolSpec.py

8/27:
- removed couple unused modules from Spec_fit_mpi.py, +..4d
- removed spectra_array=SPECTRA from Spec_fit_mpi.py, +..4d
- removed some unused code in SolSpec.py
- removed duplicate variables in SolSpec.py

9/2:
- added script using MPI to parallelize JSOC 'url' .FITS file download request

11/3:
- added shell script to download AIA (1700, 1600, 304, 171, 193, 211) and HMI (magnetogram + continuum) from JSOC using MPI
- revised JSOC download script so export request is done through shell script input
- added "-u" to shell script call of mpi_jsoc_download_complete.py, this sets the output as unbuffered, so prints correctly on Windows.  I need to check if works same on Linux.

12/21:
- added Lorentzian option to model M2 in model_sliders.py, revised script

12/26:
- added normalized averaged visual image to SolSpec.py
- added directory creation to SolSpecShell.sh (DATA directories only, not FITS)

12/27:
- added FITS directory creation to jsoc_mpi_download2.sh
- revised jsoc_mpi_download2.sh scripts to identify correct character position for file names + accept 171, 193, 211, 304, 1600, 1700, magnetogram, and continuum wavelength requests
- revised heatmap_spectra_tool.py title of plots, added white cross as current-pixel marker


2018
----

February: 
- improved heatmap spectra tool: better layout, added masking button, added button to save figures, displays masked histogram when masking heatmaps
- made heatmap spectra tool that compares different methods of weighting / uncertainties.  Displays spectra + fits in separate plot

2/22:
- added parenthesis to print statements in several scripts
- added Lorentzian model to Spec_fit_mpi.py
- revised download scripts to be compatible with Python 3

7/24:
- deleted a bunch of old / duplicate scripts
- organized the directory somewhat 