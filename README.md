# YangAndGao_AlaskaFWANT
Source data and computer codes for the Alaska FWANT paper.

## Directories
1. Figure_S2: 
Phase delay measurements for each iteration and the plotting MATLAB script.
2. Figure_S9:
Shear-wave velocities below each volcano at multiple depths.
3. Movie_S1: 
The MATLAB script to create Movie S1 in the supplement. All data, including the velocity model, volcanoes, and slab interface model, are included. A MATLAB script is provided as a wrapper to read the velocity model.
4. velocitymodel:
The 3-D shear velocity model from this study in netCDF format and a description file in plain text format. This *.nc file can be read using any netCDF reader. 

## System requirements
1. The codes have been tested on macOS HighSierra;
2. The MATLAB scripts were developed under MATLAB 2014b. However, other versions of MATLAB should also work. 
3. The MATLAB version, however, should include netCDF functions. The one used is ncread(). If not, netCDF has MATLAB package for downloads.

## Installation guide
1. Download the code to your computer: git clone https://github.com/xtyangpsp/YangAndGao_AlaskaFWANT.git
2. No further installations needed, once you have MATLAB running on your computer.

