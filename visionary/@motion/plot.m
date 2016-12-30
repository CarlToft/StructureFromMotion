function plot(m,s,color)
%MOTION/PLOT  plot 3D camera orbit (and structure)
%  plot(m,style) plots 3D camera orbit with (optional) color
%  plot(m,s,style) plots 3D camera orbit and structure with (optional) color
% INPUT:
%   m - motion
%   s - structure
%   color - plot color string

if nargin==1 | (nargin==2 & isa(s,'structure')),
  color='r';
elseif nargin==2,
  color=s;
end

holdmode=ishold;

if nargin>=2 & isa(s,'structure'),
  plot(s);
  hold on;
end


focals=focalpoints(m);
%tmp=[pflat(s),pflat(focals)];tmp(4,:)=[];ddl=std(tmp(:));
pp=pflat(getpoints(focals));

ax=axis;
ddl=std(ax)/2;

dd=zeros(3,size(m));
for i=1:size(m);
  [K,P]=rq(getcameras(m,i));
  dd(:,i)=ddl*P(3,1:3)';
end


plot(focals,[color,'*']);
hold on;axis equal;
plot(focals,[color,'o']);
plot(focals,[color,'-']);
quiver3(pp(1,:),pp(2,:),pp(3,:),dd(1,:), dd(2,:), dd(3,:),[color,'-'],'LineWidth',1.5,'MaxHeadSize',1.5);

if holdmode == 0
  hold off;
end






