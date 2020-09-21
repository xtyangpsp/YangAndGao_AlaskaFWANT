%plot 1d velocity profiles at volcanoes
%
close all;
vmodelfile='../velocitymodel/YangAndGao_AKFWANT_Vs2019_modeldata_forpaper.nc';
vtag='vs'; %velocity tag to read in the *.nc model file. 
% check the variables using ncinfo in MATLAB.
% [y,x,z,mvs]=read_netCDF_model3d(vmodelfile,vtag);
z=ncread(vmodelfile,'depth');
x=ncread(vmodelfile,'longitude');
y=ncread(vmodelfile,'latitude');
mvs=ncread(vmodelfile,vtag);
% deal with NaN values.
mvs(abs(mvs)>20)=nan;
%%
velocitytag='S';
iplot = 0; % 0, absolute velocity; 1, velocity perturbation
mvbkp=nan(size(mvs));
if strcmp(velocitytag,'P')
    mvbkp=smooth3(mvp,'box',[3 7 1]);
    %dmv=dmvp;
    figtitletag='Vp';
elseif strcmp(velocitytag,'S')
    mvbkp1=smooth3(mvs,'box',[7 15 3]);
    mvbkp2=smooth3(mvs,'box',[9 19 3]);
    mvbkp3=smooth3(mvs,'box',[11 21 3]);
    mvbkp(:,:,27:end)=mvbkp3(:,:,27:end); %below 80 km
    mvbkp(:,:,19:26)=mvbkp2(:,:,19:26); %below 40 km above 80 km
    mvbkp(:,:,1:18)=mvbkp1(:,:,1:18); %40 km to surface.
    figtitletag='Vs';
end
%
mv=mvbkp;
if iplot==1
   clear vmean3d;
   vmean3d=squeeze(mean(mean(mvbkp)));
   for i=1:length(vmean3d)
       mv(:,:,i)=100*(mvbkp(:,:,i)-vmean3d(i))/vmean3d(i);
   end
end

vmean1d=squeeze(nanmean(nanmean(mvs)));
%% get volcano locations and assocaited 1d velocity profile below it.
volcanofilebase={'ALT','DVG','BCJD','WVF'};
volcanoinfo=cell(length(volcanofilebase),1);
for i=1:length(volcanofilebase)
    filename=strcat('volcanoes_',volcanofilebase{i},'.txt');
    fidv=fopen(filename,'r');
    temp=textscan(fidv,'%f %f %*f','HeaderLines', 1);
    fclose(fidv);
    
    volcanoinfo{i}.lon=temp{1};
    volcanoinfo{i}.lat=temp{2};
    volcanoinfo{i}.tag=volcanofilebase{i};
    
    npoints=length(temp{1});
    ptval=cell(npoints,1);

    alldata=nan(length(z),npoints);
    for p=1:npoints
        ptvaltemp=nan(length(z),1);
        ptstdtemp=nan(length(z),1);
        for k=1:length(z)
            gtemp=squeeze(mv(:,:,k));
            dtemp=get_neighbor2d(x,y,gtemp,temp{1}(p),temp{2}(p),50);
            ptvaltemp(k)=nanmean(dtemp.val);
            ptstdtemp(k)=nanstd(dtemp.val);
        end
        ptval{p}.data=ptvaltemp;
        ptval{p}.std=ptstdtemp;
        ptval{p}.point=[temp{1}(p),temp{2}(p)];
        alldata(:,p)=ptvaltemp;
    end
    
    volcanoinfo{i}.vs=ptval;
    volcanoinfo{i}.vsall=alldata;    
end

%% get Moho depth
mohofile='Zhang_AlaskaHk_Moho_GRL2019_matlab_blockmean50km.txt';
mohofile4mask='Zhang_AlaskaHk_Moho_GRL2019_matlab.txt';
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

volcano_moho=cell(length(volcanofilebase),1);

for i=1:length(volcanofilebase)
    volcano_moho{i}.depth=interp2(mohox,mohoy,mohogrid,volcanoinfo{i}.lon,volcanoinfo{i}.lat);
    ptval=volcanoinfo{i}.vs;
    volcano_moho{i}.vs=nan(length(ptval),1);
    for p=1:length(ptval)
        volcano_moho{i}.vs(p)=interp1(z,ptval{p}.data,volcano_moho{i}.depth(p));
    end
end
%%

colors={'r','k',[0 0.5 0],'b'};
lines={'k-','k--','k-.','k:'};
labels={'AA','DVG','BC-JD','WVF'};
fh=[];
figure('Position',[400 400 850 800]); 
subplot(2,2,1);hold on;
for i=1:length(volcanofilebase)
    ptval=volcanoinfo{i}.vs;
    for p=1:length(ptval)
        h=plot(ptval{p}.data,z,lines{i},'color',colors{i},'linewidth',1);
    end
    fh=[fh,h];
    
end
h=plot(vmean1d,z,'k','linewidth',3.5);
fh=[fh,h];
% for i=1:length(volcanofilebase)
%     plot(volcano_moho{i}.vs,volcano_moho{i}.depth,'sw','markersize',10,...
%         'markerfacecolor',colors{i},'linewidth',2.5);
% end
hold off;
legend(fh,[labels,'Regional\newlineaverage'],'location','southwest');
ylim([15,110]);
xlabel('V_S (km/s)')
ylabel('Depth (km)')
title('(a)');
set(gca,'YDir','reverse','fontsize',14)

axis on;
box on;
grid on;

% perturbations
fh=[];
subplot(2,2,2);hold on;
for i=1:length(volcanofilebase)
    ptval=volcanoinfo{i}.vs;
    for p=1:length(ptval)
        h=plot(100*ptval{p}.data./vmean1d - 100,z,lines{i},'color',colors{i},'linewidth',1);
    end
    fh=[fh,h];
    
end
h=plot(vmean1d-vmean1d,z,'k','linewidth',3.5);
fh=[fh,h];
hold off;
legend(fh,[labels,'Regional\newlineaverage'],'location','southwest');
ylim([15,110]);
xlim([-12 12])
xlabel('dV_S relative to regional average (%)')
ylabel('Depth (km)')
title('(b)')
set(gca,'YDir','reverse','fontsize',14)

axis on;
box on;
grid on;

% the figures below plot the average across all volcanoes
fh=[];
subplot(2,2,3);hold on;
for i=1:length(volcanofilebase)
    ptval=volcanoinfo{i}.vs;
%     x1=nanmean(volcanoinfo{i}.vsall,2)+1.96*nanstd(volcanoinfo{i}.vsall,0,2);
%     x2=nanmean(volcanoinfo{i}.vsall,2)-1.96*nanstd(volcanoinfo{i}.vsall,0,2);
%     fill([reshape(x1,length(x1),1);reshape(flip(x2),length(x2),1)],[z;flip(z)],'r',...
%         'FaceColor',[0.75 .75 .75],'EdgeColor','w');
    h=plot(nanmean(volcanoinfo{i}.vsall,2),z,...
        lines{i},'color',colors{i},'linewidth',1);
    
    fh=[fh,h];
end
h=plot(vmean1d,z,'k','linewidth',3.5);
fh=[fh,h];
% for i=1:length(volcanofilebase)
%     plot(volcano_moho{i}.vs,volcano_moho{i}.depth,'sw','markersize',10,...
%         'markerfacecolor',colors{i},'linewidth',2.5);
% end
hold off;
legend(fh,[labels,'Regional\newlineaverage'],'location','southwest');
ylim([15,110]);
xlim([3.4 4.9])
xlabel('V_S (km/s)')
ylabel('Depth (km)')
title('(c)');
set(gca,'YDir','reverse','fontsize',14)

axis on;
box on;
grid on;

% perturbations
fh=[];
subplot(2,2,4);hold on;
for i=1:length(volcanofilebase)
    h=plot(100*nanmean(volcanoinfo{i}.vsall,2)./vmean1d - 100,z,lines{i},'color',colors{i},'linewidth',1);   
    fh=[fh,h];
end
h=plot(vmean1d-vmean1d,z,'k','linewidth',3.5);
fh=[fh,h];
hold off;
legend(fh,[labels,'Regional\newlineaverage'],'location','southwest');
ylim([15,110]);
xlim([-12 12])
xlabel('dV_S relative to regional average (%)')
ylabel('Depth (km)')
title('(d)')
set(gca,'YDir','reverse','fontsize',14)

axis on;
box on;
grid on;
%%
figure;
hold on;
for i=1:length(volcanofilebase)
    plot(volcanoinfo{i}.lon,volcano_moho{i}.depth,'.','color',colors{i});
end
hold off;