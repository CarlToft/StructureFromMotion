function s = structure(points,lines,quadrics)
%STRUCTURE class constructor
%  s=STRUCTURE(points,lines,quadrics) creates a STRUCTURE object where
%   points : 4xn (or 3xn) matrix where columns are homogeneous
%            (or cartesian) 3D point coordinates
%   lines  : 8xm (or 6xm) matrix where columns are two homogeneous 3D points
%            on the line (or columns are plucker coordinates)
%   quadrics : 10xp matrix where columns are dual quadric coordinates
%  All of the input arguments are optional

if nargin == 1 & isa(points,'structure'),
 s=points;
else
 if nargin >= 1,
   if size(points,1)==3
    s.points=[points;ones(1,size(points,2))];
   else
    s.points = points;
   end
 else
   s.points = zeros(4,0);
 end
 if nargin >= 2 & ~isempty(lines),
   if size(lines,1)==8,
     s.lines = lines;
   elseif size(lines,1)==6, %plucker
     zz=zeros(1,size(lines,2));
     s.lines=[zz;lines(1,:);lines(2,:);lines(3,:);...
	      -lines(1,:);zz;lines(6,:);-lines(5,:)];
   else
     error('Unknown line format');
   end
 else
   s.lines = zeros(8,0);
 end
 if nargin == 3,
   s.quadrics = quadrics;
 else
   s.quadrics = zeros(10,0);
 end

 s = class(s,'structure');

end

