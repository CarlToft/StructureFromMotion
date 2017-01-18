function xut=pextend(x);
% but=pextend(bin); add ones to get extended (or homogeneous) coordinates
% INPUT :
%  x     - matrix in which each column is a point.
%        OR structure or imagedata object with points
% OUTPUT :
%  xut    - result after normalization (in homogeneous coordinates)

if isa(x,'structure') | isa (x,'imagedata'),
  x=getpoints(x);
end

[sx,sy]=size(x);
xut=[x;ones(1,sy)];
