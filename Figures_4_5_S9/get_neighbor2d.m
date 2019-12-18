function goodpoints=get_neighbor2d(x,y,mv2d,x0,y0,maxdistxy)
%get value from neighbors in 2d grid
%By Xiaotao Yang
%Email: stcyang@gmail.com

R0=6371; %earth's radius
% %test point
% x0=-73;
% y0=44;
% z0=15;
% maxdist=15; %defined the maximum distance of the neighbors (sphere radius).
goodpoints.idx=[];
goodpoints.loc=[];
goodpoints.val=[];

% n=1;
distx=deg2km(distance(y0,x0,y,x0));
clear idx00;
idx00=find(distx<=maxdistxy)';
if ~isempty(idx00)
    for j=1:length(idx00)
        disty=deg2km(distance(y0,x0,y(idx00(j)),x));
        clear idx0;
        idx0=find(disty<=maxdistxy)';
        if ~isempty(idx0)
            goodpoints.idx=[goodpoints.idx;ones(length(idx0),1)*idx00(j),idx0];
            goodpoints.loc=[goodpoints.loc;ones(length(idx0),1)*y(idx00(j)),x(idx0)'];
            goodpoints.val=[goodpoints.val;mv2d(idx00(j),idx0)'];
        end
    end
end
    %
% if ~isempty(goodpoints.idx)
%     plot(x0,y0,'r.','markersize',25);
%     hold on;
%     plot(goodpoints.loc(:,2),goodpoints.loc(:,1),'bo','markersize',5);
%     hold off;
%     grid on;
%     box on;
% end

return;
end