function s=focalpoints(m)
%MOTION/FOCALPOINTS s=focalpoints(m) calculates camera positions
% INPUT:
%   m - motion
% OUTPUT:
%   s - structure containing 3D focal points

nbrcameras=size(m);
c=zeros(4,nbrcameras);

for i=1:nbrcameras;
  p=getcameras(m,i);
  c(:,i) = null(p);
end

s=structure(c);

