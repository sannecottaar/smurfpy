
# turfpy
toolbox utilizing receiver functions in python3

Contributors: Sanne Cottaar, Jennifer Jenkins, Stephen Pugh, Alistair Boyce, Matthew Kemp, Annemijn van Stiphout, Kieran Gilmore, and others


README last updated by: A. Boyce, 12/10/20

-----------------------------------------------------------------------
--------------------------- OUTLINE -----------------------------------
-----------------------------------------------------------------------

1. Processing_Scripts
    --> Data download, pre-processing, RF calculation, post-processing (inc. quality control)
2. Migration_Scripts
    --> Calculate RF pierce points, 1D & 3D time-to-depth conversion
3. Stacking_Scripts
    --> Epicentral distance, depth, slowness, common-conversion-point stacking
4. Plotting_Scripts
    --> Plot Pierce points, CCP stack volumes (weights, Maps, XC etc)
5. Tools
    --> Misc tools (inc velocity models)
6. Travel_Times_Slowness
    --> Reference travel time files for slowness stacking

---------------------------------------------------------------------------------
--------------------------- Processing SCRIPTS ----------------------------------
---------------------------------------------------------------------------------

# 1_download_data_per_station.py
• Description: select and download appropriate events and stations based on user inputs
• Inputs: search area lat/lon, start and end times (in function), epicentral dist of ev/station, event magnitude range, trace length, data filter band, station networks to search, dataclient
• Outputs: python stream objects in PICKLE format with a dictionary of header info for each event. Saves to ../Data/NT.STA/Originals/
• Usage: >> python3 1_download_data_per_station.py

# 2_rotate_data_NE_RT.py		
• Description: Preprocessing for RF calculation: merges truncations, trims, donsamples, rotates components Z-R-T, renames based on BAZ and EPI-DIST.
• Inputs: Data directory (usually ../Data/)
• Outputs: python stream objects in PICKLE format with a dictionary of header info for each event in a new folder leaving a copy of the unprocessed original data for future use
• Usage: >> python3 2_rotate_data_NE_RT.py

# 3_add_travel_times.py
• Description: compute predicted travel-times for user defined phases based on TauP, predicted times are added to waveform header information (python dictionary)
• Inputs: Data directory (usually ../Data/), 1D velocity model, phases to compute TTs for
• Outputs: Overwrites files from above with new dictionary (seis[0].stats.traveltimes)
• Usage: >> python3 3_add_travel_times.py P S P660s P410s

# 4_plot_data_preRF_perstation.py      
• Description: [OPTIONAL] plot Z-R-T seismogram components with predicted arrival times after initial processing, before RF creation
• Inputs: Data directory (usually ../Data/), station directory
• Outputs: On-screen plotting
• Usage: >> python3 4_plot_data_preRF_perstation.py

# receiver_function.py
• Description: Various RF algorithms: water level deconvolution, multitaper, iterative deconvolution
• Usage: called below

# 5_compute_receiver_functions.py      
• Description:
• Inputs: Data directory (usually ../Data/), horizontal component (usually radial), filter band, decon algorithm (usually iterative decon - default)
• Outputs: Adds computed RF to pre-existing PICKLE waveform file
• Usage: >> python3 5_compute_receiver_functions.py      

# 6_auto_select_receiver_functions.py  
• Description: Removes low quality ones based on set criteria: 
        1. Minimum percentage of radial compoment to be fit (after reconvolving the RF with the vertical component (fitmin)
        2. Peak amplitude max threshold before main P-wave arrival (noisebefore)
        3. Peak amplitude max threshold after main P-wave arrival (noiseafter)
        4. Peak amplitude min threshold after main P-wave arrival (minamp)
• Inputs: Data directory (usually ../Data/), horizontal component (usually radial), filter band, SNR calculation type, fitmin, noisebefore, noiseafter, minamp
• Outputs: Two ".dat" files specific to the chosen filter band recording the good RF files and the good RF file SNR ratios (V & R components)
• Usage: >> python3 6_auto_select_receiver_functions.py

# 7_plot_data_selection.py
• Description: [OPTIONAL] Plots the perstation distribution of "Acceptable - Green" and "Removed - red" events as a funciton of EQ magnitude and epicentral distance.
• Inputs: Data directory (usually ../Data/)
• Outputs: On-screen plotting
• Usage: >> python3 7_plot_data_selection.py

# 8_plot_data_perstation.py
• Description: [OPTIONAL] Plots V,R,RF as a function of time and epicentral distance.
• Inputs: Data directory (usually ../Data/), horizontal component (usually radial), filter band
• Outputs: On-screen plotting
• Usage: python3 8_plot_data_perstation.py



---------------------------------------------------------------------------------
--------------------------- MIGRATION SCRIPTS -----------------------------------
---------------------------------------------------------------------------------

# calculate_pierce_points.py
• Description: Calculate converted phase pierce points at discontinuity depths
• Inputs: Depth of piercepoints, Phase, filter band, 1D velocity model
• Outputs: Adds PP for given phase and discont depth to each Pickle file, prints to file PP_'DEPTH'km_'PHASE'_'FILTER'.txt'
• Usage: python3 calculate_pierce_points.py 410 P410s jgf1


# convert_to_depth_obspy.py
• Description: Convert RF from time to depth using 1D model (coded for Prem)
• Inputs: Filter band, 1D velocity model
• Outputs: Adds dictionary seis[0].conversions['<nameof1Dmodel>'] to each Pickle file
• Usage: python3 convert_to_depth_obspy.py

# dep_conv_AFR_AFRP20CR_AK135.py
• [OPTIONAL] As above but based on ak135 depths, converts RF from time to depth using 3D model and appropriate crustal model
• 3D model example is AFRP20 (Boyce et al., 2020 Gcubed)
• Also accounts for 3D crustal model (See Boyce et al., 2020 supplementary material) and station elevations.



---------------------------------------------------------------------------------
---------------------------STACKING SCRIPTS -----------------------------------
---------------------------------------------------------------------------------

# epicentral_distance_stack.py
• Description: Stacks the RFs in bins of epicentral distance to show the most prominent features
• Inputs: Bin size, smoothing, histogram, pierce point depth, alt/lon box, predicted travel time curves, epicentral distance range, filter band
• Outputs: Epicentral distance stack plot
• Usage: python3 epicentral_distance_stack.py


# depth_stack.py
• Description: Stacks all the RFs within the bounds for the depth stated producing one trace.
• Inputs: Data directory, filter band, conversion, normalisation depth & factor, pierce point depth and polygon.
• Outputs: Depth stack pickle file and pdf/png
• Usage: python3 depth_stack.py

# slowness_stack.py
• Description: Plot of slowness against time, using a specfic epicentral reference distance.
• Inputs: Data directory, filter band, epicentral distance nounds and reference, pierce point depth and polygon.
• Outputs: Slowness stack pickle file and pdf/png
• Usage: python3 slowness_stack.py

# depth_slowness_stack.py
• [OPTIONAL] Description: Combines previously calculated depth and slowness stack in one figure.
• Inputs: Depth and slowness stack pickle files
• Outputs: Combined depth and slowness stack image pdf/png
• Usage: python3 depth_slowness_stack.py

# common_conversion_point_stack.py
• Description: Common conversion point stacking routines (see Cottaar and Deuss, 2016). Calculates weighting factor based on Lekic et al., (2011, Science) based on distance and fresnel zone at given depth 
• Inputs: see below
• Outputs:  Common conversion point stack, weights, errors
• Usage: see below

# stack_CCP.py
• Description: wrapper for the function contained within common_conversion_point_stack.py
• Inputs: Name, filter band, conversion, smoothing factor, region, NewStack, CCP volume parameters
• Outputs: Common conversion point stack
• Usage: python3 stack_CCP.py

%%%%%%%%%%%% Parallel processing of CCP stack %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% BETA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# stack_CCP_par_beta.py
# common_conversion_point_stack_par_beta.py
• [OPTIONAL] Scripts as above but a beta version coded in Parallel.
• Computes CCP-subvolume for each RF in each station (specify max cores) and sums to master volume.

---------------------------------------------------------------------------------
--------------------------- Plotting SCRIPTS ------------------------------------
---------------------------------------------------------------------------------

# plot_map_pierce_points.py
• Description: Plots discontinuity depth pierce points
• Inputs: discontinuity depth, converted phase
• Outputs: matplotlib plot
• Usage: python3 plot_map_pierce_points.py


# CCP_plottingroutines.py
• Description: Routines for various CCP stack plot types (discontinuity depth/sampling maps, cross sections, moveout)
• Inputs: see below
• Outputs: Varoius matplotlib plot windows.
• Usage: see below

# plot_CCP.py
• Description: Wrapper for the function contained within CCP_plottingroutines.py
• Inputs: Name, filter band, conversion, smoothing factor
• Outputs: Various matplotlib plot windows.
• Usage: python3 plot_CCP.py

---------------------------------------------------------------------------------
----------------------  Travel_Times_Slowness  ----------------------------------
---------------------------------------------------------------------------------

Multiple reference travel-time files used in slowness stacking.

---------------------------------------------------------------------------------
----------------------         Tools           ----------------------------------
---------------------------------------------------------------------------------

# /MODELS/

Various models used for 1D and 3D depth stacking examples

# /PLOTTING/

# Africa_AFRP20_RF_CR1.py
• Description: Script used to plot 3D tomographic model of Boyce et al., 2020 in 3D time to depth conversion example
• Inputs: Plot type, various model parameters
• Outputs: Pdf tomogrpahic model plot
• Usage: python3 Africa_AFRP20_RF_CR1.py

---------------------------------------------------------------------------------
----------------------  South_Africa_Data      ----------------------------------
---------------------------------------------------------------------------------

Test data set for 45 stations in the XA network (doi:10.7914/SN/XA_1997)
This data is unprocessed. 
To use, copy directory to a new directory called "Data"
Then proceed from script Processing_Scripts/2_rotate_data_NE_RT.py 

---------------------------------------------------------------------------------
---------------------------------- References -----------------------------------
---------------------------------------------------------------------------------

Boyce, A. and Bastow, I. D. and Cottaar, S. and Kounoudis, R. and Guilloud De Courbeville, J. and Caunt, E. and Desai, S. (2020: Manuscript under review at Geochemistry, Geophysics, Geosystems) AFRP20: New P-wavespeed Model for the African Mantle Reveals Two Whole-Mantle Plumes Below East Africa and Neoproterozoic Modification of the Tanzania Craton

Cottaar, S., & Deuss, A.(2016).Large-scale mantle discontinuity topography beneath Europe: Signature of akimotoite in subducting slabs.J. Geophys. Res., 121(1), 279–292. doi:10.1002/2015jb012452


Lekic, V., French, S. W., & Fischer, K. M.    (2011).    Lithospheric Thinning Beneath Rifted Regions of Southern California. Science, 334 (6057), 783–787. doi:10.1126/science.1208898


