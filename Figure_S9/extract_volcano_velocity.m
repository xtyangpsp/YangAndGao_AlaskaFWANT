%analyze volcano distribution
clear vdatatmp vdata slabdata vtext;
[vdatatmp, vtext]=xlsread('AKvolclatlong.xlsx','Query1');
nl=size(vdatatmp,1);
vdata=sortrows(vdatatmp(2:nl,[4,3,8,6]),1); %sort selected columns by longitudes

% get velocity for at different depths: 80 km, 51 km, and 28 km
load 'AlaskaLayeredVelocity_griddepth_S_ite_0.05deg_05.mat';

%% map area
maparea.lon=[-158,-139]; % temporary 
maparea.lat=[57, 65];
% idxsub=find(vdata(:,1)>=maparea.lon(1) & vdata(:,1)<=maparea.lon(2));

clear idxt_ALT idxt_NEALT idxt_BCJD idxt_WVF;
idxt_ALT=find(vdata(:,1)<= -150 & vdata(:,1) >= maparea.lon(1) & ...
    vdata(:,2) >= maparea.lat(1) & vdata(:,2) <= 62.5);
% idxt_NEALT=find(vdata(:,1) <= -150 & vdata(:,1) >= -154.5 & ...
%     vdata(:,2) >= 59.8 & vdata(:,2) <= 62.5);
idxt_BCJD=find(vdata(:,1) <= -146 & vdata(:,1) >= -150 & ...
    vdata(:,2) >= 63 & vdata(:,2) <= 65);
idxt_WVF=find(vdata(:,1) <= -140 & vdata(:,1) >= -146 & ...
    vdata(:,2) >= 61 & vdata(:,2) <= 63);

%% plot Vs v.s. volcano locations (longitude)
figure('Position',[400 400 700 450]);
figlabel={'a ','b ','c ','d ','e ','f '};
clear idxnotnan;
% markersize=7;
capsize=0.5;
barlinewidth=0.25;
markersize=7;
psizefactor=1.4;

% plot locations.
subplot(2,3,1);
hold on;
plot(vdata(idxt_ALT,1),vdata(idxt_ALT,2),...
    'r^','markersize',markersize);
% plot(vdata(idxt_NEALT,1),vdata(idxt_NEALT,2),...
%     'wo','markersize',markersize,'markeredgecolor',[0 0.5 0]);
plot(vdata(idxt_BCJD,1),vdata(idxt_BCJD,2),...
    'wo','markersize',markersize,'markeredgecolor',[0 0.5 0]);
plot(vdata(idxt_WVF,1),vdata(idxt_WVF,2),...
    'bs','markersize',markersize);
hold off;
legend('AA','BC-JD','WVF','location','southeast');
xlim(maparea.lon);
ylim(maparea.lat);
box on;
axis on;
grid on;
xlabel('Longitude (degree)');
ylabel('Latitude (degree)');
title('a');
set(gca,'fontsize',14,'TickDir','out');
drawnow;

%vout is a struct, with three members: data, depth, tag
plotidx=2;
nvolvano=size(vdata,1);
% maxdist=[20,30,40,50,60];
% depths=[20 28 51 80 92];
maxdist=[50,50,50,50,50];
refmod=[3.8 3.8 4.5 4.5 4.5]; %reference model from IASP91
averagemod=[3.8 3.9 4.5 4.4 4.3]; %regional average from our model.
ylimarray=[3.2  4.4;3.2  4.3;3.8  5.1;4  4.7;4  4.7];
gxx=squeeze(XSIM(:,1,1)); %XSIM and YSIM are model grid mash for latitude and longitude, respectively.
gyy=squeeze(YSIM(1,:,1))-360;
% gzz=squeeze(ZSIM(1,1,:));
for i=2:4
    disp(['Plotting ',vout{i}.tag,' ...']);
    meanvs=nan(nvolvano,1);
    minvs=nan(nvolvano,1);
    maxvs=nan(nvolvano,1);
    
    for j=1:nvolvano
        clear dtmp;
        dtmp=get_neighbor2d(gyy,gxx,vout{i}.data,vdata(j,1),vdata(j,2),maxdist(i));
        if ~isempty(dtmp.val)
            meanvs(j)=nanmean(dtmp.val);
            minvs(j)=nanmin(dtmp.val);
            maxvs(j)=nanmax(dtmp.val);
        end
    end
    clear idxnotnan;
    idxnotnan=find(~isnan(meanvs));
    subplot(2,2,plotidx);
    hold on;
        
    %plot reference model from IASP91
    hh1=plot(maparea.lon,[refmod(i) refmod(i)],'m--','linewidth',2);
    hh2=plot(maparea.lon,[averagemod(i) averagemod(i)],'m-','linewidth',1);
    
    h1=terrorbar(vdata(intersect(idxt_ALT,idxnotnan),1),meanvs(intersect(idxt_ALT,idxnotnan)),...
        meanvs(intersect(idxt_ALT,idxnotnan))-minvs(intersect(idxt_ALT,idxnotnan)), ...
        maxvs(intersect(idxt_ALT,idxnotnan))-meanvs(intersect(idxt_ALT,idxnotnan)),capsize,'units');
    set(h1,'LineWidth',barlinewidth,'Color','r');
    plot(vdata(intersect(idxt_ALT,idxnotnan),1),meanvs(intersect(idxt_ALT,idxnotnan)),...
        'r^','markersize',markersize,'color','r','markeredgecolor','r');
%     h2=errorbar(vdata(intersect(idxt_NEALT,idxnotnan),1),meanvs(intersect(idxt_NEALT,idxnotnan)),...
%         meanvs(intersect(idxt_NEALT,idxnotnan))-minvs(intersect(idxt_NEALT,idxnotnan)), ...
%         maxvs(intersect(idxt_NEALT,idxnotnan))-meanvs(intersect(idxt_NEALT,idxnotnan)),...
%         'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    h2=terrorbar(vdata(intersect(idxt_BCJD,idxnotnan),1),meanvs(intersect(idxt_BCJD,idxnotnan)),...
        meanvs(intersect(idxt_BCJD,idxnotnan))-minvs(intersect(idxt_BCJD,idxnotnan)), ...
        maxvs(intersect(idxt_BCJD,idxnotnan))-meanvs(intersect(idxt_BCJD,idxnotnan)),capsize,'units');
    set(h2,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
    plot(vdata(intersect(idxt_BCJD,idxnotnan),1),meanvs(intersect(idxt_BCJD,idxnotnan)),...
        'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    h3=terrorbar(vdata(intersect(idxt_WVF,idxnotnan),1),meanvs(intersect(idxt_WVF,idxnotnan)),...
        meanvs(intersect(idxt_WVF,idxnotnan))-minvs(intersect(idxt_WVF,idxnotnan)), ...
        maxvs(intersect(idxt_WVF,idxnotnan))-meanvs(intersect(idxt_WVF,idxnotnan)),capsize,'units');
    set(h3,'LineWidth',barlinewidth,'Color','b');
    plot(vdata(intersect(idxt_WVF,idxnotnan),1),meanvs(intersect(idxt_WVF,idxnotnan)),...
        'bs','markersize',markersize,'markeredgecolor','b');

    legend([hh1 hh2],'IASP91','Regional');
    hold off;
    xlim(maparea.lon);
    ylim(ylimarray(i,:));
    box on;
    axis on;
    grid on;
    xlabel('Longitude (degree)');
    ylabel('Vs (km/s)');
    title([figlabel{plotidx}, ' Vs at ',vout{i}.tag]);
    set(gca,'fontsize',14,'TickDir','out');
    drawnow;
    
    plotidx=plotidx+1;
end