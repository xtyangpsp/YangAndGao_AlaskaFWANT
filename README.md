# YangAndGao_AlaskaFWANT
Source data and computer codes for the Alaska FWANT paper.

## Directories
1. Figure_S3:
Phase delay measurements distribution (histograms) for each iteration and the plotting MATLAB script.
2. Figure_6_7_S5:
Shear-wave velocities below each volcano at multiple depths. The slab interface model file is saved under Movie_S1, which can be copied or softlinked here.
3. Movie_S1:
The MATLAB script to create Movie S1 in the supplement. All data, including the velocity model, volcanoes, and slab interface model, are included. A MATLAB script is provided as a wrapper to read the velocity model under the folder: velocitymodel. The netCDF model file is also under the folder: velocitymodel.
4. velocitymodel:
The 3-D shear velocity model from this study in netCDF format and a description file in plain text format. This *.nc file can be read using any netCDF reader. The MATLAB function read_netCDF_model3d.m is used to read the model file.

## System requirements
1. The codes have been tested on macOS HighSierra;
2. The MATLAB scripts were developed under MATLAB 2014b and 2019b. However, other versions of MATLAB should also work.
3. The MATLAB version, however, should include netCDF functions. The one used is ncread(). If not, netCDF has MATLAB package for downloads.

## Installation guide
1. Download the code to your computer: git clone https://github.com/xtyangpsp/YangAndGao_AlaskaFWANT.git
2. No further installations needed, once you have MATLAB running on your computer.
