%%
% This script reads in netCDF velocity model in 3D.
%
% Prepared by Xiaotao Yang 
% Contact: stcyang@gmail.com

function [lat,lon, depth,vm]=read_netCDF_model3d(ncmodelfile,vtag)
    depth=nc_varget(ncmodelfile,'depth');
    lon=nc_varget(ncmodelfile,'longitude');
    lat=nc_varget(ncmodelfile,'latitude');
    
    vm=nc_varget(ncmodelfile,vtag);
end