function s=addlines(s,lines)
%STRUCTURE/ADDLINES s=addlines(s,lines) adds 3D lines
% Each line is respresented  by two 3D homogeneous points
% (or alternatively, by plucker coordinates)

if size(lines,1)==8,
  s.lines = [s.lines,lines];
  
elseif size(lines,1)==6, %plucker
  zz=zeros(1,size(lines,2));
  s.lines=[s.lines,[zz;lines(1,:);lines(2,:);lines(3,:);...
	      -lines(1,:);zz;lines(6,:);-lines(5,:)]];
else
  error('Unknown line format');
end

