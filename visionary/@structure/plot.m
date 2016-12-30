function plot(s,p,l)
%STRUCTURE/PLOT plot 3D structure
% plot(s,p,l)
% INPUT:
%   s - structure
%   p - point style (optional)
%   l - line style (optional)


if nargin<=1 | isempty(p),
  p = '*';
end
if nargin<=2 | isempty(l),
  l = '-';
end

holdmode=ishold;

if holdmode == 0,
  clf;
  axis([0 1 0 1 0 1]);
end

if size(s.points,2)>0,
  plot3( (s.points(1,:)./s.points(4,:))',...
	(s.points(2,:)./s.points(4,:))',...
        (s.points(3,:)./s.points(4,:))',p);
  hold on;
end
if size(s.quadrics,2)>0,
  plotquadric(s.quadrics);
  hold on;
end
if size(s.lines,2)>0,
  plotline(s.lines,l);
  hold on;
end
if holdmode == 0
  hold off;
end






