Step-1 (optional): run S1_extract_volcano_matrix.m to extract the needed data for each volcano.
Step-2: run S2_plot_volcano_matrix.m to plot the figures in the paper.

Data files needed:
* AKvolcano_geochem_namecalibrated.xlsx - geochemistry data for the volcanic rocks, downloaded from Alaska Volcano Observatory's geochemical database (https://avo.alaska.edu/geochemold).
* AKvolcanoes_extractedmatrix.xlsx - output pf Step-1 after converting to xlsx format from csv. 
* AKvolclatlong_subset.xlsx - volcano information, including locations, volumes, heat flow, and ages.
* AlaskaLayeredVelocity_griddepth_S_ite_0.05deg_05.mat - shear-wave velocity at the selected depths as shown in Figure 4 of the manuscript. This is used to provide the velocity data in extracting the velocities below each volcano at different depths.
* JFB_AKQDatabase_merged.xlsx - Heat flow data from Batir, J. F., Blackwell, D. D., Richards, M. C. Heat flow and temperature-depth curves throughout Alaska: finding regions for future geothermal exploration. J Geophys Eng, 13(3): 366-U278 (2016). Merged for different data sources. 
* S1_extract_volcano_matrix.m - MATLAB script for step-1.
* S2_plot_volcano_matrix.m - MATLAB script for step-2.
* SlabE125_ready.dat - Slab interface depth model E125 by Jadamec, M. A., Billen, M. I. Reconciling surface plate motions with rapid three-dimensional mantle flow around a slab edge. Nature, 465(7296): 338-U389 (2010).
* Zhang_AlaskaHk_Moho_GRL2019_matlab.txt - Moho depth data using teleseismic P-wave receiver functions, provided by Zhang, Y., Li, A. B., Hu, H. Crustal Structure in Alaska From Receiver Function Analysis. Geophys. Res. Lett., 46:  (2019).
* Zhang_AlaskaHk_Moho_GRL2019_matlab_blockmean50km.txt - Same as "Zhang_AlaskaHk_Moho_GRL2019_matlab.txt" but averaged with 50x50 km grids.
* volcanoes_ALT.txt - This file and the following 3 files contain volcano location data and slab depth below each volcanoes grouped by geographic areas. AA: Aleutian arc volcanoes (triangles); DVG: Denali volcanic gap (star); BC-JD: Buzzard Creek-Jumbo Dome volcanoes (circles); WVF: Wrangell volcanic field (squares). 
* volcanoes_BCJD.txt
* volcanoes_DVG.txt
* volcanoes_WVF.txt
