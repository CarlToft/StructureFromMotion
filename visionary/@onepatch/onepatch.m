function s=onepatch(points,intensities,a);
% IMCONIC/IMCONIC constructor
% s=imconic(points,normals,stddevs);
%

if nargin == 0,
  s.points = [];
  s.intensities = [];
  s.a = [];
  s = class(s,'onepatch');
elseif isa(points,'onepatch');
  s=patch;
else
  if nargin<2,
   s.points  = points;
   s.intensities = zeros(1,size(points,2));
   s.a = 1;
  else
   s.points  = points;
   s.intensities = intensities;
   s.a = a;
  end;
  s = class(s,'onepatch');
end;
