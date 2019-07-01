# YangAndGao_NatComm2019
Codes and data for Yang and Gao, Nature Communications, 2019 paper

## Directories
1. Figure_S2: 
Phase delay measurements for each iteration and the plotting MATLAB script.
2. Movie_S1: 
The MATLAB script to create Movie S1 in the supplement. All data, including the velocity model, volcanoes, and slab interface model, are included.
3. velocitymodel: 
The 3-D shear velocity model from this study in netCDF format. This *.nc file can be read using any netCDF reader. A MATLAB script is provided as a wrapper to read the velocity model.

## System requirements
1. The codes have been tested on macOS HighSierra;
2. The MATLAB scripts were developed under MATLAB 2014b. However, other versions of MATLAB should also work. 
3. The MATLAB version, however, should include netCDF functions. The one used is ncread(). If not, netCDF has MATLAB package for downloads.

## Installation guide
1. Download the code to your computer: git clone https://github.com/xtyangpsp/YangAndGao_NatComm2019.git
2. No further installations needed, once you have MATLAB running on your computer.

