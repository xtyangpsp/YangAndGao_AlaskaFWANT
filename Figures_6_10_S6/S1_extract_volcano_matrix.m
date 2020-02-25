% extract needed information (the following) and save them to excel
%     1) slab depth
%     2) volcano groups
%     3) shear velocity
%     4) major geochemistry: SiO2, MgO* (?; not implemented here yet)
%     5) volcano eruption age code
%     6) eruption volume
%     7) heatflow

%load volcano information
clear indatatable indatasubset slabdata;
%%
indatafile='AKvolclatlong_subset.xlsx';
outdatafile='AKvolcanoes_extractedmatrix.xlsx';
indatatable=readtable(indatafile,'Sheet','Query1','Range','A1:J97');
% [vdatatmp, vtext]=xlsread('AKvolclatlong_subset.xlsx','Query1');
nl=size(indatatable,1);
indatasubset=indatatable(:,[1,2,3,4,6,8,10]);
vlon=table2array(indatasubset(:,4)); %convert longitude and latitude for easy access later.
vlat=table2array(indatasubset(:,3));
slabdata=load('SlabE125_ready.dat');
% map area
maparea.lon=[-158,-139]; % 
maparea.lat=[57, 65];
% load pre-computed layer velocities:
load 'AlaskaLayeredVelocity_griddepth_S_ite_0.05deg_05.mat';
mohofile='Zhang_AlaskaHk_Moho_GRL2019_matlab_blockmean50km.txt';
mohofile4mask='Zhang_AlaskaHk_Moho_GRL2019_matlab.txt';

%heatflow
heatflowfile='JFB_AKQDatabase_merged.xlsx';
%% mask slab depth grid
clear slabgrid slabx slaby slabmask;
slabx=min(slabdata(:,1)):0.1:max(slabdata(:,1));
slaby=min(slabdata(:,2)):0.05:max(slabdata(:,2));
slabgrid=griddata(slabdata(:,1),slabdata(:,2),slabdata(:,3),slabx,slaby');
% this step is slow
didx=boundary(slabdata(:,1),slabdata(:,2),1);
doutline.lon=slabdata(didx,1);
doutline.lat=slabdata(didx,2);
%mask slab
slabmask=nan(size(slabgrid));
for i=1:size(slabmask,1)
        clear id00;
        id00=inpolygon(slabx,slaby(i)*ones(size(slabmask,2),1),...
            doutline.lon,doutline.lat);
        if ~isempty(id00);slabmask(i,id00)=1;end
end
slabgrid=slabmask.*slabgrid;

%
%
clear datat;
datat=load(mohofile);
mohox=min(datat(:,1)):1:max(datat(:,1));
mohoy=min(datat(:,2)):0.5:max(datat(:,2));
mohogrid=griddata(datat(:,1),datat(:,2),datat(:,3),mohox,mohoy');

clear moho;
moho=load(mohofile4mask);
mohoidx=boundary(moho(:,1),moho(:,2),1);
mohooutline.lon=moho(mohoidx,1);
mohooutline.lat=moho(mohoidx,2);
%mask moho
mohomask=nan(size(mohogrid));
for i=1:size(mohomask,1)
        clear id00;
        id00=inpolygon(mohox,mohoy(i)*ones(size(mohomask,2),1),...
            mohooutline.lon,mohooutline.lat);
        if ~isempty(id00);mohomask(i,id00)=1;end
end
mohogrid=mohomask.*mohogrid;

% heatflow grid
%
clear heatdata heatx heaty heatgrid;
% heatdata=load(heatflowfile);
heatdata=xlsread(heatflowfile);
heatx=min(heatdata(:,1)):.1:max(heatdata(:,1));
heaty=min(heatdata(:,2)):0.05:max(heatdata(:,2));
heatgrid=griddata(heatdata(:,1),heatdata(:,2),heatdata(:,9),heatx,heaty');

heatidx=boundary(heatdata(:,1),heatdata(:,2),1);
heatoutline.lon=heatdata(heatidx,1);
heatoutline.lat=heatdata(heatidx,2);
%mask heatflow data
heatmask=nan(size(heatgrid));
for i=1:size(heatmask,1)
        clear id00;
        id00=inpolygon(heatx,heaty(i)*ones(size(heatmask,2),1),...
            heatoutline.lon,heatoutline.lat);
        if ~isempty(id00);heatmask(i,id00)=1;end
end
heatgrid=heatmask.*heatgrid;
% get the slab depth below each volcano, through 2-D inerpolation.
slabdepth=interp2(slabx,slaby,slabgrid,vlon,vlat);

% get Moho depth below each volcano
mohodepth=interp2(mohox,mohoy,mohogrid,vlon,vlat);

% get heatflow at each volcano
heatflow=interp2(heatx,heaty,heatgrid,vlon,vlat);
indatasubset=[indatasubset table(slabdepth,mohodepth,heatflow)];
% vdataslabdepth=sortrows([vdata,vslabdepth,vmohodepth],1); %sort by longitude
%% get slab gradient
[gx,gy]=gradient(slabgrid); %gradient by default assumes the element spacing is one. 
% We will correct the gradient after extracting the values for each
% volcanoes.
R0=6371;
dx=deg2km(diff(slabx(1:2)),R0*cos(degtorad(slabx)));
dy=deg2km(diff(slaby(1:2)))*ones(size(slaby));

for i=1:size(gx,1)
    gx(i,:)=gx(i,:)./dx;
end
for i=1:size(gy,2)
    gy(:,i)=gy(:,i)./dy';
end
gxy=sqrt(gx.^2 + gy.^2); %get total gradient

slabdipgrid=radtodeg(atan(gxy)); %get dip angle at each grid point in degrees.

slabdip=interp2(slabx,slaby,slabdipgrid,vlon,vlat);
indatasubset=[indatasubset table(slabdip)];
%% set group index
% idxsub=find(vlon>=maparea.lon(1) & vlat<=maparea.lon(2));

clear idxt_ALT idxt_DVG idxt_BCJD idxt_WVF;
idxt_ALT=find(vlon<= -150 & vlon >= maparea.lon(1) & ...
    vlat >= maparea.lat(1) & vlat <= 62.5);
idxt_DVG=find(vlon <= -150.5 & vlon >= -151 & ...
    vlat >= 62.8 & vlat <= 63.2);
idxt_BCJD=find(vlon <= -146 & vlon >= -150 & ...
    vlat >= 63 & vlat <= 65);
idxt_WVF=find(vlon <= -140 & vlon >= -146 & ...
    vlat >= 61 & vlat <= 63);

group=cell(nl,1);
for i=1:length(idxt_ALT);group{idxt_ALT(i)}='AA'; end
for i=1:length(idxt_DVG);group{idxt_DVG(i)}='DVG'; end
for i=1:length(idxt_BCJD);group{idxt_BCJD(i)}='BC-JD'; end
for i=1:length(idxt_WVF);group{idxt_WVF(i)}='WVF'; end

indatasubset=[indatasubset table(group)];
%save volcanoes by areas
fid1=fopen('volcanoes_ALT.txt','w');
fprintf(fid1,'Lon Lat SlabDep\n');
for i=1:length(idxt_ALT)
    fprintf(fid1,'%g %g %g\n',vlon(idxt_ALT(i)),vlat(idxt_ALT(i)),...
        slabdepth(idxt_ALT(i)));
end
fclose(fid1);
fid2=fopen('volcanoes_DVG.txt','w');
fprintf(fid2,'Lon Lat SlabDep\n');
for i=1:length(idxt_DVG)
    fprintf(fid2,'%g %g %g\n',vlon(idxt_DVG(i)),vlat(idxt_DVG(i)),...
        slabdepth(idxt_DVG(i)));
end
fclose(fid2);
fid3=fopen('volcanoes_BCJD.txt','w');
fprintf(fid3,'Lon Lat SlabDep\n');
for i=1:length(idxt_BCJD)
    fprintf(fid3,'%g %g %g\n',vlon(idxt_BCJD(i)),vlat(idxt_BCJD(i)),...
        slabdepth(idxt_BCJD(i)));
end
fclose(fid3);
fid4=fopen('volcanoes_WVF.txt','w');
fprintf(fid4,'Lon Lat SlabDep\n');
for i=1:length(idxt_WVF)
    fprintf(fid4,'%g %g %g\n',vlon(idxt_WVF(i)),vlat(idxt_WVF(i)),...
        slabdepth(idxt_WVF(i)));
end
fclose(fid4);

%% extract velocity below volcanoes.
indatasubset_temp=indatasubset; %copy in case of errors.
clear idxnotnan;
% markersize=7;
%vout is a struct, with three members: data, depth, tag

nvolvano=nl;
% maxdist=[20,30,40,50,60];
% depths=[20 28 51 80 92];
maxdist=[50,50,50,50,50];
refmod=[3.8 3.8 4.5 4.5 4.5]; %reference model from IASP91
averagemod=[3.8 3.9 4.5 4.4 4.3]; %regional average from our model.
ylimarray=[3.2  4.4;3.2  4.3;3.8  5.1;4  4.85;4  4.7];
gxx=squeeze(XSIM(:,1,1)); %XSIM and YSIM are model grid mash for latitude and longitude, respectively.
gyy=squeeze(YSIM(1,:,1))-360;
% gzz=squeeze(ZSIM(1,1,:));
for i=2:4
    disp(['Extracting ',vout{i}.tag,' ...']);
    Vs_mean=nan(nvolvano,1);
    Vs_median=nan(nvolvano,1);
    Vs_std=nan(nvolvano,1);
    
    for j=1:nvolvano
        clear dtmp;
        dtmp=get_neighbor2d(gyy,gxx,vout{i}.data,vlon(j),vlat(j),maxdist(i));
        if ~isempty(dtmp.val)
            Vs_mean(j)=nanmean(dtmp.val);
            Vs_median(j)=nanmedian(dtmp.val);
            Vs_std(j)=nanstd(dtmp.val);
        end
    end
    
    indatasubset=[indatasubset table(Vs_mean,Vs_median,Vs_std,...
        'VariableNames',{strcat('Vs_mean_',num2str(mean(vout{i}.depth)),'km'),...
        strcat('Vs_median_',num2str(mean(vout{i}.depth)),'km'),...
        strcat('Vs_std_',num2str(mean(vout{i}.depth)),'km')})];
end

%% SiO2
%The relationships between SiO2 content and crustal melting temperature and
%fractions are from Till et al., Nature Communications, 2019
%The melt fraction is from 0 (0%) to 1 (100%), providing the estimate of
%melt fraction necessary to explain the observed composition (SiO2 content
%here).
crustalmelt_temperature = @(x) -0.5945*x.^2 + 61.054*x - 565.1;
crustalmelt_fraction = @(x) -0.0012*x.^2 + 0.1148*x - 2.1525;

%%
clear vdatachem;
vdatachem=readtable('AKvolcano_geochem_namecalibrated.xlsx','Sheet','ReadyVersion','Range','A1:P7095');
% below the row of 7095, there is no volcano name, after sorting by volcano
% names.

% subset
clear vdatachem_sub_tmp;
vdatachem_sub=vdatachem(:,{'Volcano','Longitude','Latitude','SiO2'});
for i=4
    for j=1:height(vdatachem_sub)
        if abs(table2array(vdatachem_sub(j,i)) - -9999.9)<0.001
            vdatachem_sub(j,i)=array2table(nan);
        end
    end
end

%% group by volcanoes
vchemlon=table2array(vdatachem_sub(:,2));
vchemlat=table2array(vdatachem_sub(:,3));
idxt_chemall=find(vchemlon<= -130 & vchemlon >= -169);
clear chemgroup vdatachem_sub_tmp;
vdatachem_sub_tmp=vdatachem_sub(idxt_chemall,:);

dSiO2=vdatachem_sub_tmp.SiO2;
meltT=crustalmelt_temperature(dSiO2);
meltf=crustalmelt_fraction(dSiO2);
%
vdatachem_sub_tmp=[vdatachem_sub_tmp table(meltT,meltf)];
chemgroup=sortrows(grpstats(vdatachem_sub_tmp,'Volcano',{'median','mean','std'}),'mean_Longitude','descend');

%% search for chemistry data for volcanoes in the original catalog
SiO2_mean=nan(length(indatasubset.Volcano),1);
SiO2_median=nan(length(indatasubset.Volcano),1);
SiO2_std=nan(length(indatasubset.Volcano),1);
meltT_mean=nan(length(indatasubset.Volcano),1);
meltT_median=nan(length(indatasubset.Volcano),1);
meltT_std=nan(length(indatasubset.Volcano),1);
meltf_mean=nan(length(indatasubset.Volcano),1);
meltf_median=nan(length(indatasubset.Volcano),1);
meltf_std=nan(length(indatasubset.Volcano),1);

for i=1:length(indatasubset.Volcano)
    vname=indatasubset.Volcano{i};
    idx0=find(strcmp(vname,chemgroup.Volcano)==1);
    if ~isempty(idx0)
        SiO2_std(i)=chemgroup.std_SiO2(idx0);
        SiO2_mean(i)=chemgroup.mean_SiO2(idx0);
        SiO2_median(i)=chemgroup.median_SiO2(idx0);
        
        meltT_std(i)=chemgroup.std_meltT(idx0);
        meltT_mean(i)=chemgroup.mean_meltT(idx0);
        meltT_median(i)=chemgroup.median_meltT(idx0);
        
        meltf_std(i)=chemgroup.std_meltf(idx0);
        meltf_mean(i)=chemgroup.mean_meltf(idx0);
        meltf_median(i)=chemgroup.median_meltf(idx0);
    else
        warning(['[ ',vname,' ] not found. Skipped!']);
    end
end
tabletemp=table(SiO2_std,SiO2_mean,SiO2_median,meltT_std,meltT_mean,meltT_median,meltf_std,meltf_mean,meltf_median);
%%
indatasubset=[indatasubset tabletemp];
% save table to xlsx. MAC version Excel doesn't allow saving from table. Change to csv instead.
writetable(indatasubset,strrep(outdatafile,'xlsx','csv'));