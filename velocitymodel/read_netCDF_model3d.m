%%
% This script reads in netCDF velocity model in 3D.
%
% Prepared by Xiaotao Yang 
% Contact: stcyang@gmail.com

function [lat,lon, depth,vm]=read_netCDF_model3d(ncmodelfile,vtag)
    depth=ncread(ncmodelfile,'depth');
    lon=ncread(ncmodelfile,'longitude');
    lat=ncread(ncmodelfile,'latitude');
    
    vm=ncread(ncmodelfile,vtag);
end