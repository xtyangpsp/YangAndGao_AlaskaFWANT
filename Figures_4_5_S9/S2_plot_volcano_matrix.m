% plotting volcano related analysis
%     1) GMT map of volcanoes, color-coded by slab depth. show circles around the volcanoes to illustrate the area used in extracting velocities.
%     2) Plot the velocity at multiple depths for all groups with respect to longitude.
%     3) velocity v.s. slab depth
%     4) velocity v.s. volcano age
%     5) velocity v.s. eruption volume
%     6) velocity v.s. heatflow
%     7) velocity v.s. SiO2
%% load master excel for volcano matrix
% clear all
datafile='AKvolcanoes_extractedmatrix.xlsx';
% map area
maparea.lon=[-158,-139]; % 
maparea.lat=[57, 65];
indata=readtable(datafile);
%fill empty group name with 'other'
for i=1:length(indata.group)
    if isempty(indata.group{i})
        indata.group{i}='other';
    end
end
descentrate=sin(degtorad(indata.slabdip)).*indata.PlateMotion_MORVEL;
indata=[indata table(descentrate)];
%% GMT plot of volcanoes, color-coded by slab depth
gmtmaparea=[num2str(maparea.lon(1)) '/' num2str(maparea.lon(2)) '/' ...
              num2str(maparea.lat(1)) '/' num2str(maparea.lat(2))];
gmtmapsize=[num2str(mean(maparea.lon)) '/' num2str(mean(maparea.lat)) '/' ...
    num2str(maparea.lat(1)) '/' num2str(maparea.lat(2)) '/2.8i'];
psfilenm='Figure_4a_Alaska_volcanoes_byareas.ps';
% cptfile='GMT_globe2.cpt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot using GMT matlab API to single PS file for each layer.
% To create a Grid structure from a 2-D Z array and a 1x9 header vector:
% G = gmt ('wrapgrid', Z, head)
% header is a vector with [x_min x_max, y_min y_max z_min z_max reg x_inc y_inc]
% reg: registration. we use 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear vG vcpt gheader;
% vcpt=['myjetwr_' num2str(80) '_' num2str(160) '.cpt'];
gmt('gmtset MAP_FRAME_TYPE = plain');
gmt('gmtset MAP_FRAME_WIDTH 0.1c');
gmt('gmtset FONT_TITLE = 14p,Helvetica,black');
gmt('gmtset MAP_TITLE_OFFSET = 0p');
% gmt(['makecpt -Cjetwr -T' num2str(80) '/' num2str(160) ...
%     '/',num2str(10),' > ' vcpt]);
%gmt('destroy');

gmt(['pscoast -R' gmtmaparea ' -JL' gmtmapsize ' -Di -N2/.5p,150 ',...
    '-N1/2p,80 -W.5p,80 -G200 -A50+r -Bf3a6/f1a2:."Volcanoes"' ':SEWn -X4 -Y13 -P -K  > ' psfilenm]);
%gmt('destroy');
circlesize='0.6';
gmt(['psxy -R -JL -h1 volcanoes_ALT.txt -Sc',circlesize,...
  '  -W1p,0 -O -K >> ', psfilenm]);
gmt(['psxy -R -JL -h1 volcanoes_DVG.txt -Sc',circlesize,...
  '  -W1p,0 -O -K >> ', psfilenm]);
gmt(['psxy -R -JL -h1 volcanoes_BCJD.txt -Sc',circlesize,...
  '  -W1p,0 -O -K >> ', psfilenm]);
gmt(['psxy -R -JL -h1 volcanoes_WVF.txt -Sc',circlesize,...
  '  -W1p,0 -O -K >> ', psfilenm]);

gmt(['psxy -R -JL -h1 volcanoes_ALT.txt -St0.27 ',...
  '  -W1p,red -O -K >> ', psfilenm]);
gmt(['psxy -R -JL -h1 volcanoes_DVG.txt -Sa0.3 ',...
  '  -W1p,0 -O -K >> ', psfilenm]);  
gmt(['psxy -R -JL -h1 volcanoes_BCJD.txt -Sc0.27 ',...
  '  -W1p,0/128/0 -O -K >> ', psfilenm]);    
gmt(['psxy -R -JL -h1 volcanoes_WVF.txt -Ss0.27 ',...
  '  -W1p,blue -O >> ', psfilenm]);

% COLOR BY SLAB DEPTH
% gmt(['psxy -R -JL -h1 volcanoes_ALT.txt -St0.27 ',...
%   '  -W.5p,55 -C',vcpt,' -O -K >> ', psfilenm]);
% gmt(['psxy -R -JL -h1 volcanoes_DVG.txt -Sa0.3 ',...
%   '  -W.5p,55 -C',vcpt,' -O -K >> ', psfilenm]);  
% gmt(['psxy -R -JL -h1 volcanoes_BCJD.txt -Sc0.27 ',...
%   '  -W.5p,55 -C',vcpt,' -O -K >> ', psfilenm]);    
% gmt(['psxy -R -JL -h1 volcanoes_WVF.txt -Ss0.27 ',...
%   '  -W.5p,55 -C',vcpt,' -O -K >> ', psfilenm]);
% gmt(['psscale -D3.7/-1/4.7/0.25h -C',vcpt,' -B20 -O >> ' psfilenm]);

gmt('destroy');

disp(['File saved to: ',psfilenm]);

%% get group index
clear idxt_ALT idxt_DVG idxt_BCJD idxt_WVF;
idxt_ALT=[];
idxt_DVG=[];
idxt_BCJD=[];
idxt_WVF=[];
for i=1:length(indata.group)
   if strcmp(indata.group(i),'AA')
       idxt_ALT=[idxt_ALT;i];
   elseif strcmp(indata.group(i),'DVG')
       idxt_DVG=[idxt_DVG;i];
   elseif strcmp(indata.group(i),'BC-JD')
       idxt_BCJD=[idxt_BCJD;i];
   elseif strcmp(indata.group(i),'WVF')
       idxt_WVF=[idxt_WVF;i];
   end
end

% clear statarray;
% indata_subset=indata(:,{'Longitude','Latitude','Volume_km3','EruptionCode','slabdepth',...
%     'mohodepth','heatflow','group','meanvs_28km','medianvs_28km','minvs_28km','maxvs_28km',...
%     'meanvs_51km','medianvs_51km','minvs_51km','maxvs_51km','meanvs_80km','medianvs_80km',...
%     'minvs_80km','maxvs_80km'});
% statarray=grpstats(indata_subset,'group',{'min','median','mean','max'});

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
plot(indata.Longitude(idxt_ALT),indata.Latitude(idxt_ALT),...
    'r^','markersize',markersize);
plot(indata.Longitude(idxt_DVG),indata.Latitude(idxt_DVG),...
    'wp','markersize',1.3*markersize,'markeredgecolor',[0 0 0]);
plot(indata.Longitude(idxt_BCJD),indata.Latitude(idxt_BCJD),...
    'wo','markersize',markersize,'markeredgecolor',[0 0.5 0]);
plot(indata.Longitude(idxt_WVF),indata.Latitude(idxt_WVF),...
    'bs','markersize',markersize);
hold off;
legend('AA','DVG','BC-JD','WVF','location','southeast');
xlim(maparea.lon);
ylim(maparea.lat);
daspect([1 cosd(mean(maparea.lat)) 1]);
box on;
axis on;
% grid on;
xlabel('Longitude (degree)');
ylabel('Latitude (degree)');
title('a');
set(gca,'fontsize',14,'TickDir','out');
drawnow;
%
%vout is a struct, with three members: data, depth, tag
plotidx=2;

depthlist=[80 51 28];
refmod=[4.5 4.5 3.8]; %reference model from IASP91
averagemod=[4.4 4.5 3.9]; %regional average from our model.
vlimarray=[4.1  4.8;4  4.85;3.2  4.2];

vstype='mean';
for i=1:length(depthlist)
    disp(['Plotting ',num2str(depthlist(i)),' km ...']);
    plotvs=table2array(indata(:,['Vs_' vstype '_' num2str(depthlist(i)) 'km']));
%     medianvs=table2array(indata(:,['medianvs_' num2str(depthlist(i)) 'km']));
    stdvs=table2array(indata(:,['Vs_std_' num2str(depthlist(i)) 'km']));
%     maxvs=table2array(indata(:,['maxvs_' num2str(depthlist(i)) 'km']));
    
    clear idxnotnan;
    idxnotnan=find(~isnan(plotvs));
    subplot(2,2,plotidx);
    hold on;
        
    %plot reference model from IASP91
    
    hh1=plot(maparea.lon,[averagemod(i) averagemod(i)],'k-','color',[0.5 .5 .5],'linewidth',2);
    hh2=plot(maparea.lon,[refmod(i) refmod(i)],'k--','color',[0.2 .2 .2],'linewidth',2.5);
    
    h1=terrorbar(indata.Longitude(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
        stdvs(intersect(idxt_ALT,idxnotnan)),stdvs(intersect(idxt_ALT,idxnotnan)),capsize,'units');
    set(h1,'LineWidth',barlinewidth,'Color','r');
    h1p=plot(indata.Longitude(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
        'r^','markersize',markersize,'color','r','markeredgecolor','r');
    h2=terrorbar(indata.Longitude(intersect(idxt_DVG,idxnotnan)),plotvs(intersect(idxt_DVG,idxnotnan)),...
        stdvs(intersect(idxt_DVG,idxnotnan)), stdvs(intersect(idxt_DVG,idxnotnan)),capsize,'units');
        set(h2,'LineWidth',barlinewidth,'Color',[0 0 0]);
    h2p=plot(indata.Longitude(intersect(idxt_DVG,idxnotnan)),plotvs(intersect(idxt_DVG,idxnotnan)),...
            'kp','markersize',1.3*markersize,'color',[0 0 0],'markeredgecolor',[0 0 0]);
    h3=terrorbar(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
        stdvs(intersect(idxt_BCJD,idxnotnan)), stdvs(intersect(idxt_BCJD,idxnotnan)),capsize,'units');
    set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
    h3p=plot(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
        'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    h4=terrorbar(indata.Longitude(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
        stdvs(intersect(idxt_WVF,idxnotnan)),stdvs(intersect(idxt_WVF,idxnotnan)),capsize,'units');
    set(h4,'LineWidth',barlinewidth,'Color','b');
    h4p=plot(indata.Longitude(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
        'bs','markersize',markersize,'markeredgecolor','b');

    if i==3
        legend([hh1 hh2],'Regional average','IASP91','location','southwest');
    elseif i==2
        legend([h1p h2p h3p h4p],'AA','DVG','BC-JD','WVF','location','northwest');
    end
    text(maparea.lon(2)-7,vlimarray(i,2)-0.1,['depth: ',num2str(depthlist(i)), ' km'],'fontsize',14)
    hold off;
    xlim(maparea.lon);
    ylim(vlimarray(i,:));
    box on;
    axis on;
%     grid on;
    xlabel('Longitude (degree)');
    ylabel('Vs (km/s)');
    title([figlabel{plotidx}]);
    set(gca,'fontsize',14,'TickDir','out');
    drawnow;
    
    plotidx=plotidx+1;
end
figname=['Figure_4_Alaska_volcanoes_',vstype,'vs.eps'];
% set(gcf,'PaperPositionMode','auto');
saveas(gca,figname,'epsc');
disp(['Saved to: ' figname]);

%% plot Vs perturbation v.s. volcano locations (longitude)
figure('Position',[400 400 850 200]);
figlabel={'a ','b ','c ','d ','e ','f '};
clear idxnotnan;
% markersize=7;
capsize=0.5;
barlinewidth=0.25;
markersize=7;
psizefactor=1.4;

%
%vout is a struct, with three members: data, depth, tag
plotidx=1;

depthlist=[80 51 28];
refmod=[4.5 4.5 3.8]; %reference model from IASP91
averagemod=[4.4 4.5 3.9]; %regional average from our model.
dvlimarray=[-7  10;-12  7;-17  6];

vstype='mean';
for i=1:length(depthlist)
    disp(['Plotting ',num2str(depthlist(i)),' km ...']);
    plotvs=table2array(indata(:,['Vs_' vstype '_' num2str(depthlist(i)) 'km']));
    plotvs=100*(plotvs - averagemod(i))/averagemod(i);
%     medianvs=table2array(indata(:,['medianvs_' num2str(depthlist(i)) 'km']));
    stdvs=table2array(indata(:,['Vs_std_' num2str(depthlist(i)) 'km']));
    stdvs=100*stdvs/averagemod(i);
%     maxvs=table2array(indata(:,['maxvs_' num2str(depthlist(i)) 'km']));
    
    clear idxnotnan;
    idxnotnan=find(~isnan(plotvs));
    subplot(1,3,plotidx);
    hold on;
        
    %plot
    plot(maparea.lon,[0 0],'k');
    h1=terrorbar(indata.Longitude(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
        stdvs(intersect(idxt_ALT,idxnotnan)),stdvs(intersect(idxt_ALT,idxnotnan)),capsize,'units');
    set(h1,'LineWidth',barlinewidth,'Color','r');
    h1p=plot(indata.Longitude(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
        'r^','markersize',markersize,'color','r','markeredgecolor','r');
    h2=terrorbar(indata.Longitude(intersect(idxt_DVG,idxnotnan)),plotvs(intersect(idxt_DVG,idxnotnan)),...
        stdvs(intersect(idxt_DVG,idxnotnan)), stdvs(intersect(idxt_DVG,idxnotnan)),capsize,'units');
        set(h2,'LineWidth',barlinewidth,'Color',[0 0 0]);
    h2p=plot(indata.Longitude(intersect(idxt_DVG,idxnotnan)),plotvs(intersect(idxt_DVG,idxnotnan)),...
            'kp','markersize',1.3*markersize,'color',[0 0 0],'markeredgecolor',[0 0 0]);
    h3=terrorbar(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
        stdvs(intersect(idxt_BCJD,idxnotnan)), stdvs(intersect(idxt_BCJD,idxnotnan)),capsize,'units');
    set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
    h3p=plot(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
        'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    h4=terrorbar(indata.Longitude(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
        stdvs(intersect(idxt_WVF,idxnotnan)),stdvs(intersect(idxt_WVF,idxnotnan)),capsize,'units');
    set(h4,'LineWidth',barlinewidth,'Color','b');
    h4p=plot(indata.Longitude(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
        'bs','markersize',markersize,'markeredgecolor','b');

    if i==2
        legend([h1p h2p h3p h4p],'AA','DVG','BC-JD','WVF','location','south');
    end
    text(maparea.lon(2)-8,dvlimarray(i,2)-1.5,['depth: ',num2str(depthlist(i)), ' km'],'fontsize',12)
    
    if i==3
        text(maparea.lon(1)+1,dvlimarray(i,1)+3.5,['average'],'fontsize',12)
        text(maparea.lon(1)+1,dvlimarray(i,1)+1.5,[num2str(averagemod(i)),' km/s'],'fontsize',12)
    else
        text(maparea.lon(1)+1,dvlimarray(i,2)-1.5,['average'],'fontsize',12)
        text(maparea.lon(1)+1,dvlimarray(i,2)-3,[num2str(averagemod(i)),' km/s'],'fontsize',12)
    end
    hold off;
    xlim(maparea.lon);
    ylim(dvlimarray(i,:));
    box on;
    axis on;
%     grid on;
    xlabel('Longitude (degree)');
    ylabel('Vs anomaly (%)');
    title([figlabel{plotidx}]);
    set(gca,'fontsize',12,'TickDir','out');
    drawnow;
    
    plotidx=plotidx+1;
end
figname=['Figure_S9_Alaska_volcanoes_',vstype,'dvs.eps'];
% set(gcf,'PaperPositionMode','auto');
saveas(gca,figname,'epsc');
disp(['Saved to: ' figname]);
%% plots of velocity v.s. [slab depth, eruption volume, volcano age, heatflow, SiO2] at three depths
fontsize=12;
figure('Position',[1400 200 800 850]);
figlabel={'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r'};
capsize=0.4;
barunits='units';
barlinewidth=0.25;
markersize=7;
psizefactor=1.4;
depthlist=[80 51 28];
averagemod=[4.4 4.5 3.9]; %regional average from our model at the corresponding depths.
varlist={'slabdepth','Volume_km3','EruptionCode','heatflow','SiO2_mean'}; %'heatflow','PlateMotion_MORVEL'
varlabellist={'Slab depth (km)','log_1_0 (Volume in km^3)',...
    'Last eruption age code','Heat flow (mW/m^2)','SiO_2 (%)'}; %'Heat flow (mW/m^2)','Plate motion (mm/year)'
titlelist={'Slab depth','Volume','Last eruption age code','Heat flow','SiO_2'}; %'Heat flow','Plate motion'
vstype='mean';
vlimarray=[4.1  4.8;4  4.85;3.3  4.2];
% xlimarray=[0 3.1;0  9;60 180;92 97.5];
% textposx=[1.8, 5.3, 130, 95.2];
xlimarray=[60 180; 0 3.1;0 9;92 97.5;47 73];%;50 65
textposx=[65, .2,.6,92.3,61];%, 57
subplotpar=[length(varlist) length(depthlist)];
%
fit_depth=nan(length(varlist)*length(depthlist),1);
fit_par=cell(length(varlist)*length(depthlist),1);
fit_p1_ALT=nan(length(varlist)*length(depthlist),1);
fit_p0_ALT=nan(length(varlist)*length(depthlist),1);
fit_R2_ALT=nan(length(varlist)*length(depthlist),1);
fit_p1_BCJD=nan(length(varlist)*length(depthlist),1);
fit_p0_BCJD=nan(length(varlist)*length(depthlist),1);
fit_R2_BCJD=nan(length(varlist)*length(depthlist),1);
fit_p1_WVF=nan(length(varlist)*length(depthlist),1);
fit_p0_WVF=nan(length(varlist)*length(depthlist),1);
fit_R2_WVF=nan(length(varlist)*length(depthlist),1);
count=1;
for i=1:length(depthlist)
    disp(['Plotting ',num2str(depthlist(i)),' km ...']);
    plotvs=table2array(indata(:,['Vs_' vstype '_' num2str(depthlist(i)) 'km']));
%     medianvs=table2array(indata(:,['medianvs_' num2str(depthlist(i)) 'km']));
    stdvs=table2array(indata(:,['Vs_std_' num2str(depthlist(i)) 'km']));
    
    clear idxnotnan;
    idxnotnan0=find(~isnan(plotvs));
    for j=1:length(varlist)
        hsp=subplot(subplotpar(1),subplotpar(2),i+subplotpar(2)*(j-1));
        capsize=0.04*range(xlimarray(j,:));
        hold on;
        clear dtemp;
        dtemp=table2array(indata(:,varlist{j}));
        idxnotnan=intersect(idxnotnan0,find(~isnan(dtemp)));
        if strcmp(varlist{j},'Volume_km3'); dtemp=log10(dtemp);end
        %all four areas
%         idxall=intersect(idxnotnan,[idxt_ALT;idxt_DVG;idxt_BCJD;idxt_WVF]);
%         dtempall=sortrows([dtemp(idxall) plotvs(idxall) minvs(idxall) maxvs(idxall)],1);
%         area(dtempall(:,1),dtempall(:,4),'facecolor',[0.7 .7 .7],'edgecolor','none');
%         area(dtempall(:,1),dtempall(:,3),'facecolor','w','edgecolor','none');
%         plot(dtempall(:,1),dtempall(:,2),'k-','color',[0.6 .6 .6],'linewidth',5);
        if strcmp(varlist{j},'SiO2_mean')
            area([53,65,65,53],[vlimarray(i,1),vlimarray(i,1),...
                vlimarray(i,2),vlimarray(i,2)],'facecolor',[0.7 .7 .7],'edgecolor','none')
%             plot([53 53],[vlimarray(i,1) vlimarray(i,2)],'k--');
%             plot([65 65],[vlimarray(i,1) vlimarray(i,2)],'k--');
        end
        plot(xlimarray(j,:),[averagemod(i) averagemod(i)],'k--','color',[0.2 .2 .2],'linewidth',1);
        
        p=[nan nan];
        R2=nan;
        if ~isempty(intersect(idxt_ALT,idxnotnan))
            
            h1=terrorbar(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
                stdvs(intersect(idxt_ALT,idxnotnan)),stdvs(intersect(idxt_ALT,idxnotnan)),capsize,barunits);
            set(h1,'LineWidth',barlinewidth,'Color','r');
            h1p=plot(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
                'r^','markersize',markersize,'color','r','markeredgecolor','r');
            [p,S,mu]=polyfit(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),1);
            f=polyval(p,dtemp(intersect(idxt_ALT,idxnotnan)),[],mu);
            R2=1 - (S.normr/norm(plotvs(intersect(idxt_ALT,idxnotnan)) - mean(plotvs(intersect(idxt_ALT,idxnotnan)))))^2;
            plot(dtemp(intersect(idxt_ALT,idxnotnan)),f,'k-','linewidth',1.5,'color',[0.5 0 0]);
        end
%         fprintf('%g %s %s %g %g\n',depthlist(i),varlist{j},'AA',p(1),p(2));
        fit_depth(count)=depthlist(i);
        fit_par{count}=titlelist{j};
        fit_p1_ALT(count)=p(1);
        fit_p0_ALT(count)=p(2);
        fit_R2_ALT(count)=R2;
        
        if ~isempty(intersect(idxt_DVG,idxnotnan))
            h2=terrorbar(dtemp(intersect(idxt_DVG,idxnotnan)),plotvs(intersect(idxt_DVG,idxnotnan)),...
                stdvs(intersect(idxt_DVG,idxnotnan)),stdvs(intersect(idxt_DVG,idxnotnan)),capsize,barunits);
                set(h2,'LineWidth',barlinewidth,'Color',[0 0 0]);
            h2p=plot(dtemp(intersect(idxt_DVG,idxnotnan)),plotvs(intersect(idxt_DVG,idxnotnan)),...
                    'kp','markersize',1.3*markersize,'color',[0 0 0],'markeredgecolor',[0 0 0]);
        end
        
        p=[nan nan];
        R2=nan;
        if ~isempty(intersect(idxt_BCJD,idxnotnan))
            h3=terrorbar(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
                stdvs(intersect(idxt_BCJD,idxnotnan)),stdvs(intersect(idxt_BCJD,idxnotnan)),capsize,barunits);
            set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
            h3p=plot(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
                'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
            [p,S,mu]=polyfit(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),1);
            f=polyval(p,dtemp(intersect(idxt_BCJD,idxnotnan)),[],mu);
            R2=1 - (S.normr/norm(plotvs(intersect(idxt_BCJD,idxnotnan)) - mean(plotvs(intersect(idxt_BCJD,idxnotnan)))))^2;
            plot(dtemp(intersect(idxt_BCJD,idxnotnan)),f,'g-','linewidth',1.5,'color',[0 0.4 0]);
            
        end
%         fprintf('%g %s %s %g %g\n',depthlist(i),varlist{j},'BC-JD',p(1),p(2));
        fit_p1_BCJD(count)=p(1);
        fit_p0_BCJD(count)=p(2);
        fit_R2_BCJD(count)=R2;
        
        p=[nan nan];
        R2=nan;
        if ~isempty(intersect(idxt_WVF,idxnotnan))
            h4=terrorbar(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
                stdvs(intersect(idxt_WVF,idxnotnan)),stdvs(intersect(idxt_WVF,idxnotnan)),capsize,barunits);
            set(h4,'LineWidth',barlinewidth,'Color','b');
            h4p=plot(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
                'bs','markersize',markersize,'markeredgecolor','b');
            if strcmp(varlist{j},'Volume_km3')
                dtx=dtemp(intersect(idxt_WVF,idxnotnan));
                dty=plotvs(intersect(idxt_WVF,idxnotnan));
                [p,S,mu]=polyfit(dtx(dtx>1),dty(dtx>1),1);
                f=polyval(p,dtx(dtx>1),[],mu);
                R2=1 - (S.normr/norm(dty - mean(dty)))^2;
                plot(dtx(dtx>1),f,'c-','linewidth',1.5,'color',[0 0 .5]);
            else
                [p,S,mu]=polyfit(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),1);
                f=polyval(p,dtemp(intersect(idxt_WVF,idxnotnan)),[],mu);
                R2=1 - (S.normr/norm(plotvs(intersect(idxt_WVF,idxnotnan)) - mean(plotvs(intersect(idxt_WVF,idxnotnan)))))^2;
                plot(dtemp(intersect(idxt_WVF,idxnotnan)),f,'c-','linewidth',1.5,'color',[0 0 .5]);
            end
            
        end
%         fprintf('%g %s %s %g %g\n',depthlist(i),varlist{j},'WVF',p(1),p(2));
        fit_p1_WVF(count)=p(1);
        fit_p0_WVF(count)=p(2);
        fit_R2_WVF(count)=R2;
        
        count=count+1;
        if j==length(varlist) && i==length(depthlist)
            legend([h1p h2p h3p h4p],'AA','DVG','BC-JD','WVF','location','southwest');
        end
        if strcmp(varlist{j},'SiO2_mean') && depthlist(i)==80
            text(textposx(j),vlimarray(i,1)+0.07,['depth: ',num2str(depthlist(i)), ' km'],'fontsize',12);
        else
            text(textposx(j),vlimarray(i,2)-0.07,['depth: ',num2str(depthlist(i)), ' km'],'fontsize',12);
        end
        if strcmp(varlist{j},'SiO2_mean') && depthlist(i)==80
            text(47.5,vlimarray(i,2)-0.1,'mafic','fontsize',12);
            text(54,vlimarray(i,2)-0.1,'intermediate','fontsize',12);
            text(67,vlimarray(i,2)-0.1,'silicic','fontsize',12);
        end
        hold off;
        if strcmp(varlist{j},'EruptionCode'); set(gca,'XTick',1:8);end
        xlim(xlimarray(j,:));
        ylim(vlimarray(i,:));
%         if strcmp(varlist{j},'slabdepth')
%             if depthlist(i)==80
%                 ylim([4.1 4.65])
%             elseif depthlist(i)==28
%                 ylim([3.4 4.2]);
%             end
%         end
        box on;
        axis on;
%         grid on;
        xlabel(varlabellist{j});
        ylabel('Vs (km/s)')
%         if i==1; ;end
        title([figlabel{i+subplotpar(2)*(j-1)}],'fontsize',fontsize);
        set(gca,'fontsize',fontsize,'TickDir','out');
%         a1=gca;
%         a2=axes('YAxisLocation', 'Right');
%         yt=a1.YTick;
%         set(a2, 'color', 'none')
%         set(a2, 'XTick', [])
%         a2.YTick=a1.YTick;
%         a2.YTickLabel=num2cell(-100*(1 - yt/averagemod(i)));
        drawnow;
    end
end
figname=['Figure_5_Alaska_volcanoes_',vstype,'vs_others.eps'];
% set(gcf,'PaperPositionMode','auto');
saveas(gca,figname,'epsc');
disp(['Saved to: ' figname]);

table_fit=table(fit_depth,fit_par,fit_p1_ALT,fit_p0_ALT,fit_R2_ALT,...
    fit_p1_BCJD,fit_p0_BCJD,fit_R2_BCJD,fit_p1_WVF,fit_p0_WVF,fit_R2_WVF);
%%
tablefilename=['Table_S2_Alaska_volcanoes_',vstype,'vs_others_polyfitcoeff.csv'];

writetable(table_fit,tablefilename);
%% plots of velocity v.s. the temperature and fraction of Crustal melt
figure('Position',[1400 200 800 550]);
figlabel={'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'};
capsize=0.4;
barunits='units';
barlinewidth=0.25;
markersize=7;
psizefactor=1.4;
% depthlist=[28 51 80];
% averagemod=[3.9 4.5 4.4]; %regional average from our model at the corresponding depths.

vstype='mean';
% vlimarray=[3.2  4.2;4  4.85;4.1  4.8];
xlimarray=[47 73;760 1020;10  60];
textposx=[62, 770, 12];

SiO2=table2array(indata(:,'SiO2_mean'));
SiO2_std=table2array(indata(:,'SiO2_std'));

% From Till et al., Nature Communications, 2019
%T_melt = -0.5945*SiO2^2 + 61.054*SiO2 - 565.1 %Crustal melt temperature
% f_melt = -0.0012*SiO2^2 + 0.1148*SiO2 - 2.1525 %Crustal melt fraction.

meltT=indata.meltT_mean;
meltf=indata.meltf_mean;
meltf=meltf*100;
meltT(SiO2<53)=nan;
meltf(SiO2<53)=nan;

meltT_std=indata.meltT_std;
meltf_std=indata.meltf_std;
meltf_std=meltf_std*100;
meltT_std(SiO2<53)=nan;
meltf_std(SiO2<53)=nan;

% remove data when computed crustal melt fraction - std(melt fraction) < 0
minmeltf=meltf-meltf_std;
meltT(minmeltf<0)=nan;
meltf(minmeltf<0)=nan;
SiO2(minmeltf<0)=nan;

%
subplotpar=[3 3];
for i=1:length(depthlist)
    disp(['Plotting ',num2str(depthlist(i)),' km ...']);
    plotvs=table2array(indata(:,['Vs_' vstype '_' num2str(depthlist(i)) 'km']));
%     medianvs=table2array(indata(:,['medianvs_' num2str(depthlist(i)) 'km']));
    stdvs=table2array(indata(:,['Vs_std_' num2str(depthlist(i)) 'km']));
    
    clear idxnotnan;
    idxnotnan0=find(~isnan(plotvs));
    
    % plot SiO2:
    dtemp=SiO2;
    idxnotnan=intersect(idxnotnan0,find(~isnan(dtemp)));
    subplot(subplotpar(1),subplotpar(2),i); hold on;
    capsize=0.04*range(xlimarray(1,:));
    plot([53 53],[vlimarray(i,1) vlimarray(i,2)],'k--');
    plot([65 65],[vlimarray(i,1) vlimarray(i,2)],'k--');
       
    plot(xlimarray(1,:),[averagemod(i) averagemod(i)],'k-','color',[0.5 .5 .5],'linewidth',2);
    if ~isempty(intersect(idxt_ALT,idxnotnan))
        h1=terrorbar(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
            stdvs(intersect(idxt_ALT,idxnotnan)),stdvs(intersect(idxt_ALT,idxnotnan)),capsize,barunits);
        set(h1,'LineWidth',barlinewidth,'Color','r');
        h1p=plot(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
            'r^','markersize',markersize,'color','r','markeredgecolor','r');
    end
    if ~isempty(intersect(idxt_BCJD,idxnotnan))
        h3=terrorbar(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
            stdvs(intersect(idxt_BCJD,idxnotnan)),stdvs(intersect(idxt_BCJD,idxnotnan)),capsize,barunits);
        set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
        h3p=plot(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
            'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    end
    if ~isempty(intersect(idxt_WVF,idxnotnan))
        h4=terrorbar(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
            stdvs(intersect(idxt_WVF,idxnotnan)),stdvs(intersect(idxt_WVF,idxnotnan)),capsize,barunits);
        set(h4,'LineWidth',barlinewidth,'Color','b');
        h4p=plot(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
            'bs','markersize',markersize,'markeredgecolor','b');
    end
    if i==1
        text(47.5,vlimarray(i,1)+0.1,'mafic','fontsize',12);
        text(54,vlimarray(i,1)+0.1,'intermediate','fontsize',12);
        text(67,vlimarray(i,1)+0.1,'silicic','fontsize',12);
    end
    text(textposx(1),vlimarray(i,2)-0.1,['depth: ',num2str(depthlist(i)), ' km'],'fontsize',12);
    hold off;
    xlim(xlimarray(1,:));
    ylim(vlimarray(i,:));
    box on;
    axis on;
%         grid on;
    xlabel('SiO_2 (%)');
    ylabel('Vs (km/s)');
    title([figlabel{i}]);
    set(gca,'fontsize',14,'TickDir','out');
    drawnow;
    
    %plot melt temperature
    dtemp=meltT;
    idxnotnan=intersect(idxnotnan0,find(~isnan(dtemp)));
    subplot(subplotpar(1),subplotpar(2),i+subplotpar(2)); hold on;
    capsize=0.04*range(xlimarray(2,:));
    plot(xlimarray(2,:),[averagemod(i) averagemod(i)],'k-','color',[0.5 .5 .5],'linewidth',3);
    if ~isempty(intersect(idxt_ALT,idxnotnan))
        h1=terrorbar(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
            stdvs(intersect(idxt_ALT,idxnotnan)),stdvs(intersect(idxt_ALT,idxnotnan)),capsize,barunits);
        set(h1,'LineWidth',barlinewidth,'Color','r');
        h1p=plot(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
            'r^','markersize',markersize,'color','r','markeredgecolor','r');
    end
    if ~isempty(intersect(idxt_BCJD,idxnotnan))
        h3=terrorbar(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
            stdvs(intersect(idxt_BCJD,idxnotnan)),stdvs(intersect(idxt_BCJD,idxnotnan)),capsize,barunits);
        set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
        h3p=plot(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
            'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    end
    if ~isempty(intersect(idxt_WVF,idxnotnan))
        h4=terrorbar(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
            stdvs(intersect(idxt_WVF,idxnotnan)),stdvs(intersect(idxt_WVF,idxnotnan)),capsize,barunits);
        set(h4,'LineWidth',barlinewidth,'Color','b');
        h4p=plot(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
            'bs','markersize',markersize,'markeredgecolor','b');
    end
    text(textposx(2),vlimarray(i,2)-0.1,['depth: ',num2str(depthlist(i)), ' km'],'fontsize',12);
    hold off;
    xlim(xlimarray(2,:));
    ylim(vlimarray(i,:));
    box on;
    axis on;
%         grid on;
    xlabel(['Crustal melt temperature (',char(0176),'C)']);
    ylabel('Vs (km/s)');
    title([figlabel{i+subplotpar(2)}]);
    set(gca,'fontsize',14,'TickDir','out');
    drawnow;
    
    dtemp=meltf;
    idxnotnan=intersect(idxnotnan0,find(~isnan(dtemp)));
    subplot(subplotpar(1),subplotpar(2),i+2*subplotpar(2)); hold on;
    capsize=0.04*range(xlimarray(3,:));
    plot(xlimarray(3,:),[averagemod(i) averagemod(i)],'k-','color',[0.5 .5 .5],'linewidth',3);
    if ~isempty(intersect(idxt_ALT,idxnotnan))
        h1=terrorbar(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
            stdvs(intersect(idxt_ALT,idxnotnan)),stdvs(intersect(idxt_ALT,idxnotnan)),capsize,barunits);
        set(h1,'LineWidth',barlinewidth,'Color','r');
        h1p=plot(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
            'r^','markersize',markersize,'color','r','markeredgecolor','r');
    end
    if ~isempty(intersect(idxt_BCJD,idxnotnan))
        h3=terrorbar(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
            stdvs(intersect(idxt_BCJD,idxnotnan)),stdvs(intersect(idxt_BCJD,idxnotnan)),capsize,barunits);
        set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
        h3p=plot(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
            'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    end
    if ~isempty(intersect(idxt_WVF,idxnotnan))
        h4=terrorbar(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
            stdvs(intersect(idxt_WVF,idxnotnan)),stdvs(intersect(idxt_WVF,idxnotnan)),capsize,barunits);
        set(h4,'LineWidth',barlinewidth,'Color','b');
        h4p=plot(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
            'bs','markersize',markersize,'markeredgecolor','b');
    end
    text(textposx(3),vlimarray(i,2)-0.1,['depth: ',num2str(depthlist(i)), ' km'],'fontsize',12);
    hold off;
    if i==length(depthlist)
        legend([h1p h3p h4p],'AA','BC-JD','WVF','location','northeast');
    end
    xlim(xlimarray(3,:));
    ylim(vlimarray(i,:));
    box on;
    axis on;
%         grid on;
    xlabel('Crustal melt fraction (%)');
    ylabel('Vs (km/s)');
    title([figlabel{i+2*subplotpar(2)}]);
    set(gca,'fontsize',14,'TickDir','out');
    drawnow;
end
figname=['Alaska_volcanoes_',vstype,'vs_melt.eps'];
% set(gcf,'PaperPositionMode','auto');
saveas(gca,figname,'epsc');
disp(['Saved to: ' figname]);

%% heat flow and melt fraction
heatflow=table2array(indata(:,'heatflow'));
clear idxnotnan;
idxnotnan0=find(~isnan(heatflow));
xlimarray=[45 72;780 1100;0  75];
dtemp=SiO2;
dtemp_std=SiO2_std;

idxnotnan=intersect(idxnotnan0,find(~isnan(dtemp)));
figure('Position',[1400 900 850 220]);
subplot(1,3,1); hold on;

% SiO2 along longitude
capsize=0.04*range(maparea.lon);
plot([maparea.lon(1) maparea.lon(2)],[53 53],'k--');
plot([maparea.lon(1) maparea.lon(2)],[65 65],'k--');
if ~isempty(intersect(idxt_ALT,idxnotnan))
    h1=terrorbar(indata.Longitude(intersect(idxt_ALT,idxnotnan)),dtemp(intersect(idxt_ALT,idxnotnan)),...
        dtemp_std(intersect(idxt_ALT,idxnotnan)),dtemp_std(intersect(idxt_ALT,idxnotnan)),capsize,barunits);
    set(h1,'LineWidth',barlinewidth,'Color','r');
    h1p=plot(indata.Longitude(intersect(idxt_ALT,idxnotnan)),dtemp(intersect(idxt_ALT,idxnotnan)),...
        'r^','markersize',markersize,'color','r','markeredgecolor','r');
end
if ~isempty(intersect(idxt_BCJD,idxnotnan))
    h3=terrorbar(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),dtemp(intersect(idxt_BCJD,idxnotnan)),...
        dtemp_std(intersect(idxt_BCJD,idxnotnan)),dtemp_std(intersect(idxt_BCJD,idxnotnan)),capsize,barunits);
    set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
    h3p=plot(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),dtemp(intersect(idxt_BCJD,idxnotnan)),...
        'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
end
if ~isempty(intersect(idxt_WVF,idxnotnan))
    h4=terrorbar(indata.Longitude(intersect(idxt_WVF,idxnotnan)),dtemp(intersect(idxt_WVF,idxnotnan)),...
        dtemp_std(intersect(idxt_WVF,idxnotnan)),dtemp_std(intersect(idxt_WVF,idxnotnan)),capsize,barunits);
    set(h4,'LineWidth',barlinewidth,'Color','b');
    h4p=plot(indata.Longitude(intersect(idxt_WVF,idxnotnan)),dtemp(intersect(idxt_WVF,idxnotnan)),...
        'bs','markersize',markersize,'markeredgecolor','b');
end
text(maparea.lon(2)-5,47,'mafic','fontsize',12);
text(maparea.lon(2)-8,58,'intermediate','fontsize',12);
text(maparea.lon(2)-5,70,'silicic','fontsize',12);
hold off;
xlim([maparea.lon]);
ylim(xlimarray(1,:));
% ylim([0 100])
box on;
axis on;
%         grid on;
ylabel(['SiO_2 (%)']);
xlabel('Longitude (degree)');
title([figlabel{1}]);
set(gca,'fontsize',14,'TickDir','out');
drawnow;

% Heatflow v.s. Crustal melt temperature
dtemp=meltT;
dtemp_std=meltT_std;
idxnotnan=intersect(idxnotnan0,find(~isnan(dtemp)));
subplot(1,3,2); hold on;
capsize=0.04*range(maparea.lon);
if ~isempty(intersect(idxt_ALT,idxnotnan))
    h1=terrorbar(indata.Longitude(intersect(idxt_ALT,idxnotnan)),dtemp(intersect(idxt_ALT,idxnotnan)),...
        dtemp_std(intersect(idxt_ALT,idxnotnan)),dtemp_std(intersect(idxt_ALT,idxnotnan)),capsize,barunits);
    set(h1,'LineWidth',barlinewidth,'Color','r');
    h1p=plot(indata.Longitude(intersect(idxt_ALT,idxnotnan)),dtemp(intersect(idxt_ALT,idxnotnan)),...
        'r^','markersize',markersize,'color','r','markeredgecolor','r');
end
if ~isempty(intersect(idxt_BCJD,idxnotnan))
    h3=terrorbar(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),dtemp(intersect(idxt_BCJD,idxnotnan)),...
        dtemp_std(intersect(idxt_BCJD,idxnotnan)),dtemp_std(intersect(idxt_BCJD,idxnotnan)),capsize,barunits);
    set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
    h3p=plot(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),dtemp(intersect(idxt_BCJD,idxnotnan)),...
        'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
end
if ~isempty(intersect(idxt_WVF,idxnotnan))
    h4=terrorbar(indata.Longitude(intersect(idxt_WVF,idxnotnan)),dtemp(intersect(idxt_WVF,idxnotnan)),...
        dtemp_std(intersect(idxt_WVF,idxnotnan)),dtemp_std(intersect(idxt_WVF,idxnotnan)),capsize,barunits);
    set(h4,'LineWidth',barlinewidth,'Color','b');
    h4p=plot(indata.Longitude(intersect(idxt_WVF,idxnotnan)),dtemp(intersect(idxt_WVF,idxnotnan)),...
        'bs','markersize',markersize,'markeredgecolor','b');
end

hold off;
xlim([maparea.lon]);
ylim(xlimarray(2,:));
box on;
axis on;
%         grid on;
ylabel(['Crustal melt temperature (',char(0176),'C)']);
xlabel('Longitude (degree)');
title([figlabel{2}]);
set(gca,'fontsize',14,'TickDir','out');
drawnow;

% Heatflow v.s. Crustal melt fraction (%)
dtemp=meltf;
dtemp_std=meltf_std;
idxnotnan=intersect(idxnotnan0,find(~isnan(dtemp)));
subplot(1,3,3); hold on;
capsize=0.04*range(maparea.lon);
if ~isempty(intersect(idxt_ALT,idxnotnan))
    h1=terrorbar(indata.Longitude(intersect(idxt_ALT,idxnotnan)),dtemp(intersect(idxt_ALT,idxnotnan)),...
        dtemp_std(intersect(idxt_ALT,idxnotnan)),dtemp_std(intersect(idxt_ALT,idxnotnan)),capsize,barunits);
    set(h1,'LineWidth',barlinewidth,'Color','r');
    h1p=plot(indata.Longitude(intersect(idxt_ALT,idxnotnan)),dtemp(intersect(idxt_ALT,idxnotnan)),...
        'r^','markersize',markersize,'color','r','markeredgecolor','r');
end
if ~isempty(intersect(idxt_BCJD,idxnotnan))
    h3=terrorbar(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),dtemp(intersect(idxt_BCJD,idxnotnan)),...
        dtemp_std(intersect(idxt_BCJD,idxnotnan)),dtemp_std(intersect(idxt_BCJD,idxnotnan)),capsize,barunits);
    set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
    h3p=plot(indata.Longitude(intersect(idxt_BCJD,idxnotnan)),dtemp(intersect(idxt_BCJD,idxnotnan)),...
        'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
end
if ~isempty(intersect(idxt_WVF,idxnotnan))
    h4=terrorbar(indata.Longitude(intersect(idxt_WVF,idxnotnan)),dtemp(intersect(idxt_WVF,idxnotnan)),...
        dtemp_std(intersect(idxt_WVF,idxnotnan)),dtemp_std(intersect(idxt_WVF,idxnotnan)),capsize,barunits);
    set(h4,'LineWidth',barlinewidth,'Color','b');
    h4p=plot(indata.Longitude(intersect(idxt_WVF,idxnotnan)),dtemp(intersect(idxt_WVF,idxnotnan)),...
        'bs','markersize',markersize,'markeredgecolor','b');
end

hold off;
legend([h1p,h3p,h4p],'AA','BC-JD','WVF','location','northeast');
xlim([maparea.lon]);
ylim(xlimarray(3,:));
box on;
axis on;
%         grid on;
ylabel('Crustal melt fraction (%)');
xlabel('Longitude (degree)');
title([figlabel{3}]);
set(gca,'fontsize',14,'TickDir','out');
drawnow;

figname=['Alaska_volcanoes_SiO2_melt.eps'];
% set(gcf,'PaperPositionMode','auto');
% saveas(gca,figname,'epsc');
% disp(['Saved to: ' figname]);

%%
figure('Position',[1400 200 800 220]);
idxnotnan0=find(~isnan(indata.slabdepth));
idxnotnan=intersect(idxnotnan0,find(~isnan(indata.slabdip)));
subplot(1,3,1);axis on; box on;hold on;

h1p=plot(indata.slabdip(intersect(idxt_ALT,idxnotnan)),indata.slabdepth(intersect(idxt_ALT,idxnotnan)),...
    'r^','markersize',markersize,'color','r','markeredgecolor','r');
h2p=plot(indata.slabdip(intersect(idxt_DVG,idxnotnan)),indata.slabdepth(intersect(idxt_DVG,idxnotnan)),...
        'kp','markersize',1.3*markersize,'color',[0 0 0],'markeredgecolor',[0 0 0]);
h3p=plot(indata.slabdip(intersect(idxt_BCJD,idxnotnan)),indata.slabdepth(intersect(idxt_BCJD,idxnotnan)),...
    'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
h4p=plot(indata.slabdip(intersect(idxt_WVF,idxnotnan)),indata.slabdepth(intersect(idxt_WVF,idxnotnan)),...
    'bs','markersize',markersize,'markeredgecolor','b');

hold off;
legend([h1p h2p h3p h4p],'AA','DVG','BC-JD','WVF','location','northwest');
title('a')
xlabel('Slab dip (degree)')
ylabel('Slab depth (km)')
set(gca,'fontsize',fontsize,'TickDir','out');

idxnotnan=intersect(idxnotnan0,find(~isnan(indata.PlateMotion_MORVEL)));
subplot(1,3,2);axis on; box on;hold on;

h1p=plot(indata.PlateMotion_MORVEL(intersect(idxt_ALT,idxnotnan)),indata.slabdepth(intersect(idxt_ALT,idxnotnan)),...
    'r^','markersize',markersize,'color','r','markeredgecolor','r');
h2p=plot(indata.PlateMotion_MORVEL(intersect(idxt_DVG,idxnotnan)),indata.slabdepth(intersect(idxt_DVG,idxnotnan)),...
        'kp','markersize',1.3*markersize,'color',[0 0 0],'markeredgecolor',[0 0 0]);
h3p=plot(indata.PlateMotion_MORVEL(intersect(idxt_BCJD,idxnotnan)),indata.slabdepth(intersect(idxt_BCJD,idxnotnan)),...
    'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
h4p=plot(indata.PlateMotion_MORVEL(intersect(idxt_WVF,idxnotnan)),indata.slabdepth(intersect(idxt_WVF,idxnotnan)),...
    'bs','markersize',markersize,'markeredgecolor','b');

hold off;
title('b')
xlabel('Plate motion (mm/year)')
ylabel('Slab depth (km)')
set(gca,'fontsize',fontsize,'TickDir','out');

idxnotnan=intersect(idxnotnan0,find(~isnan(descentrate)));
subplot(1,3,3);axis on; box on;hold on;

h1p=plot(descentrate(intersect(idxt_ALT,idxnotnan)),indata.slabdepth(intersect(idxt_ALT,idxnotnan)),...
    'r^','markersize',markersize,'color','r','markeredgecolor','r');
h2p=plot(descentrate(intersect(idxt_DVG,idxnotnan)),indata.slabdepth(intersect(idxt_DVG,idxnotnan)),...
        'kp','markersize',1.3*markersize,'color',[0 0 0],'markeredgecolor',[0 0 0]);
h3p=plot(descentrate(intersect(idxt_BCJD,idxnotnan)),indata.slabdepth(intersect(idxt_BCJD,idxnotnan)),...
    'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
h4p=plot(descentrate(intersect(idxt_WVF,idxnotnan)),indata.slabdepth(intersect(idxt_WVF,idxnotnan)),...
    'bs','markersize',markersize,'markeredgecolor','b');

hold off;
title('c')
xlabel('Decent rate (mm/year)')
ylabel('Slab depth (km)')
set(gca,'fontsize',fontsize,'TickDir','out');

figname=['Alaska_volcanoes_slabdepth_others.eps'];
% set(gcf,'PaperPositionMode','auto');
saveas(gca,figname,'epsc');
disp(['Saved to: ' figname]);

%% plots of velocity v.s. [slab dip, plate motion,descent rate] at three depths
fontsize=12;
figure('Position',[1400 200 800 470]);
figlabel={'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r'};
capsize=0.4;
barunits='units';
barlinewidth=0.25;
markersize=7;
psizefactor=1.4;
depthlist=[80 51 28];
averagemod=[4.4 4.5 3.9]; %regional average from our model at the corresponding depths.
varlist={'slabdip','PlateMotion_MORVEL','descentrate'}; %'heatflow','PlateMotion_MORVEL'
varlabellist={'Slab dip (degree)','Plate motion rate (mm/year)','Slab descent rate (mm/year)'}; %'Heat flow (mW/m^2)','Plate motion (mm/year)'
titlelist={'Slab dip','Plate motion rate','Slab descent rate'}; %'Heat flow','Plate motion'
vstype='mean';
vlimarray=[4.1  4.8;4  4.85;3.3  4.2];
% xlimarray=[0 3.1;0  9;60 180;92 97.5];
% textposx=[1.8, 5.3, 130, 95.2];
xlimarray=[22 55; 52 63;24 46];%;50 65
textposx=[40, 58,36];%, 57
subplotpar=[length(varlist) length(depthlist)];
%

for i=1:length(depthlist)
    disp(['Plotting ',num2str(depthlist(i)),' km ...']);
    plotvs=table2array(indata(:,['Vs_' vstype '_' num2str(depthlist(i)) 'km']));
%     medianvs=table2array(indata(:,['medianvs_' num2str(depthlist(i)) 'km']));
    stdvs=table2array(indata(:,['Vs_std_' num2str(depthlist(i)) 'km']));
    
    clear idxnotnan;
    idxnotnan0=find(~isnan(plotvs));
    for j=1:length(varlist)
        hsp=subplot(subplotpar(1),subplotpar(2),i+subplotpar(2)*(j-1));
        capsize=0.04*range(xlimarray(j,:));
        hold on;
        clear dtemp;
        
        
        if strcmp(varlist{j},'descentrate')
            dtemp=descentrate;
        else
            dtemp=table2array(indata(:,varlist{j}));
        end
        idxnotnan=intersect(idxnotnan0,find(~isnan(dtemp)));

        plot(xlimarray(j,:),[averagemod(i) averagemod(i)],'k-','color',[0.5 .5 .5],'linewidth',2);
        if ~isempty(intersect(idxt_ALT,idxnotnan))
            h1=terrorbar(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
                stdvs(intersect(idxt_ALT,idxnotnan)),stdvs(intersect(idxt_ALT,idxnotnan)),capsize,barunits);
            set(h1,'LineWidth',barlinewidth,'Color','r');
            h1p=plot(dtemp(intersect(idxt_ALT,idxnotnan)),plotvs(intersect(idxt_ALT,idxnotnan)),...
                'r^','markersize',markersize,'color','r','markeredgecolor','r');
        end
        if ~isempty(intersect(idxt_DVG,idxnotnan))
            h2=terrorbar(dtemp(intersect(idxt_DVG,idxnotnan)),plotvs(intersect(idxt_DVG,idxnotnan)),...
                stdvs(intersect(idxt_DVG,idxnotnan)),stdvs(intersect(idxt_DVG,idxnotnan)),capsize,barunits);
                set(h2,'LineWidth',barlinewidth,'Color',[0 0 0]);
            h2p=plot(dtemp(intersect(idxt_DVG,idxnotnan)),plotvs(intersect(idxt_DVG,idxnotnan)),...
                    'kp','markersize',1.3*markersize,'color',[0 0 0],'markeredgecolor',[0 0 0]);
        end
        if ~isempty(intersect(idxt_BCJD,idxnotnan))
            h3=terrorbar(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
                stdvs(intersect(idxt_BCJD,idxnotnan)),stdvs(intersect(idxt_BCJD,idxnotnan)),capsize,barunits);
            set(h3,'LineWidth',barlinewidth,'Color',[0 0.5 0]);
            h3p=plot(dtemp(intersect(idxt_BCJD,idxnotnan)),plotvs(intersect(idxt_BCJD,idxnotnan)),...
                'ko','markersize',markersize,'color',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
        end
        if ~isempty(intersect(idxt_WVF,idxnotnan))
            h4=terrorbar(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
                stdvs(intersect(idxt_WVF,idxnotnan)),stdvs(intersect(idxt_WVF,idxnotnan)),capsize,barunits);
            set(h4,'LineWidth',barlinewidth,'Color','b');
            h4p=plot(dtemp(intersect(idxt_WVF,idxnotnan)),plotvs(intersect(idxt_WVF,idxnotnan)),...
                'bs','markersize',markersize,'markeredgecolor','b');
        end
        if j==length(varlist) && i==length(depthlist)
            legend([h1p h2p h3p h4p],'AA','DVG','BC-JD','WVF','location','southwest');
        end
        text(textposx(j),vlimarray(i,2)-0.07,['depth: ',num2str(depthlist(i)), ' km'],'fontsize',12);

        hold off;

        xlim(xlimarray(j,:));
        ylim(vlimarray(i,:));

        box on;
        axis on;
%         grid on;
        xlabel(varlabellist{j});
        ylabel('Vs (km/s)')
%         if i==1; ;end
        title([figlabel{i+subplotpar(2)*(j-1)}],'fontsize',fontsize);
        set(gca,'fontsize',fontsize,'TickDir','out');

        drawnow;
    end
end
figname=['Alaska_volcanoes_',vstype,'vs_slab.eps'];
% set(gcf,'PaperPositionMode','auto');
saveas(gca,figname,'epsc');
disp(['Saved to: ' figname]);