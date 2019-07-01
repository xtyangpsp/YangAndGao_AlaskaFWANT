%%
% This script plots the 3D isosurface of the velocity model, with specified netCDF file name.
% It also makes a movie through rotating the camera position and saves it
% to mp4 format.
%
% Prepared by Xiaotao Yang 
% Contact: stcyang@gmail.com
%
% Reference:
% Yang, X., H. Gao (in revision). Seismic imaging of slab segmentation and 
% correlation with volcano distribution along the Aleutian-Alaska subduction zone, 
% Nature Communications
%
% Alaska volcano catalog was downloaded from Alaska volcano observatory (https://www.avo.alaska.edu)
%
% Slab interface model E125 was from Jadamec, M. A., and M. I. Billen (2010),
% Reconciling surface plate motions with rapid three-dimensional mantle flow 
% around a slab edge, Nature, 465(7296), 338-U389, doi:10.1038/nature09053. 

%plot isosurface
close all;

% load master data
load AlaskaBorder;
state=[]; state(1).polygon(:,1)=Alaska.lon;state(1).polygon(:,2)=Alaska.lat;

%
vmodelfile='YangAndGao_AKFWANT_Vs2019_modeldata_forpaper.nc';
vtag='vs'; %velocity tag to read in the *.nc model file. 
% check the variables using ncinfo in MATLAB.
% [y,x,z,mvs]=read_netCDF_model3d(vmodelfile,vtag);
z=ncread(vmodelfile,'depth');
x=ncread(vmodelfile,'longitude');
y=ncread(vmodelfile,'latitude');
mvs=ncread(vmodelfile,vtag);
% deal with NaN values.
mvs(abs(mvs)>20)=nan;

%
% maparea.lon=[-161,-136.5]; % temporary 
% maparea.lat=[56, 66.8];
AKvolcanoes=load('AKvolclatlong_ready_matlab.txt');
%%
clear slabdata slabgrid;
slabdata=load('SlabE125_ready.dat');
slablon=min(slabdata(:,1)):0.1:max(slabdata(:,1));
slablat=min(slabdata(:,2)):0.1:max(slabdata(:,2));
[slabX,slabY]=meshgrid(slablon,slablat);
slabgrid=griddata(slabdata(:,1),slabdata(:,2),slabdata(:,3),...
    slablon,(slablat)');
%%
vplot_all=smooth3(permute(mvs,[2,3,1]),'box',[13 25 1]); %

%% create mask grid based on seismic ray coverage for period of 25-50 s (ray path > 10)
amask=nan(length(y),length(x));
load('AlaskaRayCoverOutline_ite_0.05deg_05_25-50s_cutoff10.mat');
for i=1:size(amask,1)
    clear id00;
    id00=inpolygon(x,y(i)*ones(size(amask,2),1),raycover.data(:,1),...
            raycover.data(:,2));
    amask(i,id00)=1;
end
%
amask3d=nan(size(vplot_all));
for k=1:length(z)
    amask3d(:,:,k)=amask;
end
vplot=vplot_all.*amask3d;
%%
maparea.lon=[-158,-139]; % temporary 
maparea.lat=[57, 65];
myzlim=[30 120];
lonidx=find(x>=maparea.lon(1)+0.1 & x<=maparea.lon(2)-0.1);
latidx=find(y>=maparea.lat(1)+0.1 & y<=maparea.lat(2)-0.1);
depthidx=find(z>=myzlim(1) & z<=myzlim(2));
%
clear isov isocap_low;
zlimindexmin=min(depthidx); zlimindexmax=max(depthidx);
isovalue1=4.1;
zlimindexmax_slow=zlimindexmax;
isov=isosurface(x(lonidx),y(latidx),z(zlimindexmin+2:zlimindexmax_slow),...
    vplot(latidx,lonidx,zlimindexmin+2:zlimindexmax_slow),isovalue1);
isocap_low=isocaps(x(lonidx),y(latidx),z(zlimindexmin+2:zlimindexmax_slow),...
    vplot(latidx,lonidx,zlimindexmin+2:zlimindexmax_slow),isovalue1,'below');

isovalue=4.6;
isovUM=isosurface(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax),vplot(latidx,lonidx,zlimindexmin:zlimindexmax),isovalue);
isocap_slab=isocaps(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax),...
    vplot(latidx,lonidx,zlimindexmin:zlimindexmax),isovalue);

%% plot 4.1 km/s isosurface.
bgcolor='k'; %'w' for print; 'k' for screen play
figure('position',[1400 400 1280 720],'color',bgcolor);
surfacealpha=1;
p=patch(isov,'FaceColor',[1 0.75 0],...
   'EdgeColor','none','FaceAlpha',1);
isonormals(x,y,z,...
    vplot,p);
patch(isocap_low,'FaceColor','interp','EdgeColor','none','FaceAlpha',1)

view(3); axis tight

hold on;
myzlimUM=myzlim;%[0 150];

pUM=patch(isovUM,'FaceColor',[0 0.25 .6],...
   'EdgeColor','none','FaceAlpha',surfacealpha);

isonormals(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax),...
    vplot(latidx,lonidx,zlimindexmin:zlimindexmax),pUM);
patch(isocap_slab,'FaceColor','interp','EdgeColor','none','FaceAlpha',surfacealpha)
colormap(jetwr);
caxis([3.7 4.9]);

view(-31,25); %axis tight
camlight right
camlight left
lightangle(125,65);
lightangle(90,65); 
lighting gouraud
set(gca,'ZDir','reverse');
hold on;

% plot slab surface.
hss=surface(slabX,slabY,slabgrid,'Facecolor',[0.5 0.5 0.5],'Edgecolor',[0.5 0.5 0.5],...
    'FaceAlpha',0.8,'linewidth',.5,'EdgeAlpha',1);
if strcmp(bgcolor,'k')
    for sb=1:length(state)
        plot3(state(sb).polygon(:,1), state(sb).polygon(:,2),myzlimUM(1)*ones(length(state(sb).polygon(:,1)),1),...
            'color',[.9 .9 .9],'LineWidth',2);
    end
else
    for sb=1:length(state)
        plot3(state(sb).polygon(:,1), state(sb).polygon(:,2),myzlimUM(1)*ones(length(state(sb).polygon(:,1)),1),...
            'color',[.2 .2 .2],'LineWidth',2);
    end
end

plot3(AKvolcanoes(:,1),AKvolcanoes(:,2),myzlimUM(1)*ones(length(AKvolcanoes(:,2)),1),...
    '^','linewidth',1,'markersize',10,'markeredgecolor',[.9 .9 .9],'markerfacecolor',[1 0 1]);

axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2) myzlimUM(1) myzlimUM(2)]);

hold off

box on;
ax = gca;
ax.BoxStyle = 'full';
% grid on;
daspect([1 cosd(mean(maparea.lat)) 12]);
%
hcbar=colorbar('eastoutside');
set(hcbar,'TickDirection','out','Ticks',3.7:.2:4.9);

hcbar.Label.String='V_S (km/s)';

if strcmp(bgcolor,'k')
    xlabel('Longitude','fontsize',16,'color','w');
    ylabel('Latitude','fontsize',16,'color','w');
    zlabel('Depth (km)','fontsize',16,'color','w');
    set(gca,'Fontsize',16,'ZTick',myzlimUM(1):30:myzlimUM(2),'ZTickLabel',myzlimUM(1):30:myzlimUM(2),...
        'YTick',maparea.lat(1):2:maparea.lat(2),'XTick',round(maparea.lon(1):3:round(maparea.lon(2)))); 
    ax.XColor='w';
    ax.YColor='w';
    ax.ZColor='w';
    hcbar.Color='w';
    hcbar.Label.Color='w';
    set(gca,'Color','k')
else
    xlabel('Longitude','fontsize',16);
    ylabel('Latitude','fontsize',16);
    zlabel('Depth (km)','fontsize',16);
    set(gca,'Fontsize',16,'ZTick',myzlimUM(1):30:myzlimUM(2),'ZTickLabel',myzlimUM(1):30:myzlimUM(2),...
        'YTick',maparea.lat(1):2:maparea.lat(2),'XTick',round(maparea.lon(1):3:round(maparea.lon(2))));
%zlim([40 100]);
end
%% making movie.
savevideo=1;
if savevideo
    aviid = VideoWriter([velocitytag,'_isosurface_',...
        num2str(isovalue1),'_and_',num2str(isovalue),'_',bgcolor,'_withcolorbar.mp4'],'MPEG-4');
    open(aviid);
end
view(-31,25);

pause(0.5)
for i=1:10
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for el=25:.5:41
    view(-31,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for az=-31:-.5:-59
    view(az,41);
    camzoom(1.005);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
view(az,43);
pause(1);
for i=1:25
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end

for az=-59:.5:-9
    el=el+.25;
    view(az,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end

for az=-9:.5:14
%     el=el+.5;
    view(az,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
pause(1);
for i=1:25
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for cz=1:20
    camzoom(1.019);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
pause(1);
for i=1:25
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for az=14:1:85
%     el=el+.5;
    view(az,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
pause(1);
for i=1:25
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for cz=1:20
    camzoom(1.012);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
pause(1);
for i=1:25
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
%
for az=85:1:120
    el=el+.2;
    view(az,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
%
pause(1);
for i=1:25
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for cz=1:20
    camzoom(1/1.016);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
%
for i=1:15
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
% pause;
for az=120:1:160
    el=el-.2;
    view(az,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
pause(1);
for i=1:25
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
%
for az=160:1:180
    el=el+.3;
    view(az,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
% pause
for az=-180:1:-130
    el=el-.4;
    view(az,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for i=1:20
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for az=-130:1:-60
    el=el-.4;
    view(az,el);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
pause(1);
for i=1:25
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
%
for az=-60:.35:-25
    el=el+.15;
    view(az,el);
    camzoom(1/1.004);
    pause(.08);
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
for i=1:20
    if savevideo
        F = getframe(gcf);
        writeVideo(aviid,F);
    end
end
if savevideo;close(aviid);end