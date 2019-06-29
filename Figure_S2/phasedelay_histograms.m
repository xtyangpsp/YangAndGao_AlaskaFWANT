%%
% This script plots the phase delay histograms as shwon in Figure S2 in
% the supplementary information.
%
% Prepared by Xiaotao Yang 
% Contact: stcyang@gmail.com
%
% Reference:
% Yang, X., H. Gao (in revision). Seismic imaging of slab segmentation and 
% correlation with volcano distribution along the Aleutian-Alaska subduction zone, 
% Nature Communications

close all; clear all;

tempdelaydatafilelist={'ite_0.05deg_01_measure_result_Alaska_simpleaverage.dat_latlon',...
    'ite_0.05deg_02_measure_result_Alaska.dat_latlon',...
    'ite_0.05deg_03_measure_result_Alaska.dat_latlon',...
    'ite_0.05deg_04_measure_result_Alaska.dat_latlon',...
    'ite_0.05deg_05_measure_result_Alaska.dat_latlon',...
    'ite_0.05deg_06_measure_result_Alaska.dat_latlon'};
% tempdelaydatafilelist={'ite03_measure_result_NEANT_finalclean.dat_latlon'};
numinput=numel(tempdelaydatafilelist); 
%%%% define parameters %%%%
max_dT = 10; % maximum delay in second (see mesure_phase_delay.csh)

%Alaska fband
pband=[100, 200; 75,150;50,100;35,75;25,50;15,35;10,25;7.5,15];fband=flip(1./pband,2);
[nfb, nc]=size(fband);
      
%% 
figure('Position',[400 400 1050 650]);
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) ','(i) '};
colorlist={'b','r','r','g','r','m'};
cmap=jet(2*numinput);
dtbinwidth=0.5;
dtbin=-max_dT-0.5*dtbinwidth:dtbinwidth:max_dT+0.5*dtbinwidth;
minnum=1000000*ones(nfb,1);maxnum=zeros(nfb,1);
stdall=nan(numinput,nfb);
for k=1:numinput
    tempdelaydatafile=tempdelaydatafilelist{k};
    fidtemp=fopen(tempdelaydatafile,'r');
    tempdelaydata=textscan(fidtemp,'%s %s %f %f %f %f %f %f %s %f %*f');
    fclose(fidtemp);
    [src,rcv,slat,slon,rlat,rlon,delay,err,fb,mycoe] = tempdelaydata{1:10};
    clear tempdelaydata;
    nd = length(slat);

    % get values for each frequence band
    idf=nan(length(delay),nfb); %idf: index of frequency band, initiated as full length;
    dataf.delay=idf;
    dataf.coe=idf;
    % dataf.rms=nan(nfb,1);
    dataf.num=nan(nfb,1);
    % dataf.src=cell(length(delay),nfb);
    % dataf.rcv=cell(length(delay),nfb);

    for i=1:nfb
        disp(strcat(num2str(pband(i,1)),'-',num2str(pband(i,2)),' s'));
        ftag=strcat('f',num2str(i));
        idftemp=strmatch(ftag,fb);
        dataf.num(i)=length(idftemp);
        idf(1:dataf.num(i),i)=idftemp;

        dataf.delay(1:dataf.num(i),i)=delay(idftemp);
        dataf.coe(1:dataf.num(i),i)=mycoe(idftemp);
    %     dataf.rms(i)=rms(dataf.delay(1:dataf.num(i),i));

    %     for j=1:dataf.num(i)
    %         dataf.src{j,i}=src(idftemp(j));
    %         dataf.rcv{j,i}=rcv(idftemp(j));
    %     end
    end

    %
    disp('Plotting ...');

    for i=1:nfb
        %plot
        if dataf.num(i) >0
            clear idx;
            idx=find(abs(dataf.delay(:,i)) <= max_dT & dataf.coe(:,i) >= 0.75);
%             length(idx);
            stdall(k,i)=std(dataf.delay(idx,i));
            meanall(k,i)=mean(dataf.delay(idx,i));
            subplot(3,3,i), hold on, box on
            
            [n,xout]=histcounts(dataf.delay(idx,i),dtbin);
            if min(n)< minnum(i); minnum(i)=min(n);end
            if max(n)> maxnum(i); maxnum(i)=max(n);end
            
%             bar(xout, n, 'EdgeColor',colorlist{k},'FaceColor','none','LineWidth',.5);
            plot(dtbin(1:end-1)+0.5*dtbinwidth, n, '-','linewidth',2,'color',...
                cmap(2*k-1,:));
            xlim([-max_dT-.5 max_dT+.5]);
%             text(-max_dT+1,max(n)*0.8,['\sigma=' num2str(stdall(k,i),2)],'color',cmap((k-1)*55+1,:));
%             c = length(xout);
%             w = xout(c)-xout(1);
%             t = linspace(xout(1)-w/c,xout(end)+w/c,c+1);
%             dt = diff(t);
%             Fvals = cumsum([0,n.*dt]);
%             F = spline(t, [0, Fvals, 0]);
%             DF = fnder(F);  % computes its first derivative
%     %         fnplt(DF, colorlist{k}, 2);
% 
%             maxfnval(k,i)=max(fnval(DF,xout(1):.1:xout(end)));
            if k==numinput
                plot([0 0],[minnum(i) 1.1*maxnum(i)],'k-','LineWidth',1)
                ylim([minnum(i) 1.1*maxnum(i)]);
                %text(max_dT-1,100,['RMS = ' num2str(dataf.rms(i),3)]);
                xlabel('Delay (s)','FontSize',14);
                ylabel('Count','FontSize',14);
                set(gca,'XTick',-max_dT:2:max_dT,'XTickLabel',-max_dT:2:max_dT);
                if i >3 && i <6
                   xlim([-5 5]);
                    set(gca,'XTick',-5:1:5,'XTickLabel',-5:1:5); 
                elseif i >=6
                   xlim([-5 5]);
                    set(gca,'XTick',-5:1:5,'XTickLabel',-5:1:5); 
                end
                set(gca,'Fontsize',14,'TickDir','out');
    %             title(strcat(figlabel{i}, num2str(int16(pband(i,1))),'-',num2str(int16(pband(i,2))),' s:  ',...
    %                 num2str(dataf.num(i))),'FontSize',14);
                title([figlabel{i},num2str(pband(i,1)),'-',...
                    num2str(pband(i,2)),' s'],'FontSize',14);
            end
            %pause;
        end
    end
    
    drawnow;
end 
% subplot(2,4,1);legend('0.025-1','0.025-2','0.025-3','0.015-1','0.015-2','0.015-3','0.015-4');
% subplot(3,3,1);legend('reference','preferred');
subplot(3,3,1);legend('ref.','ite-1','ite-2','ite-3','ite-4','ite-5','location','northwest');


subplot(3,3,9);
% plot(stdall(1,:),'o-','color',cmap((1-1)*55+1,:),'linewidth',2*2.1^(1-1)); hold on;
% plot(stdall(2,:),'^-','color',cmap((2-1)*55+1,:),'linewidth',2*2.1^(2-1));
errorbar(meanall(1,:),stdall(1,:),'o','color',cmap(1,:),'linewidth',1); hold on;
errorbar(meanall(numinput,:),stdall(numinput,:),'.','color',cmap(2*numinput-1,:),'linewidth',2,'markersize',20);
ylabel('mean (s)','FontSize',14);
xlim([.5 nfb+.5]);
ylim([-8.5 8.5]);
plot([.5 nfb+.5],[0 0],'k-','LineWidth',1);
set(gca,'XTick',1:1:nfb,'XTickLabel',figlabel,'TickDir','out');
legend('reference','preferred');
set(gca,'Fontsize',14,'YTick',-8:4:8);
title([figlabel{9} 'average delays'],'FontSize',14);
hold off;