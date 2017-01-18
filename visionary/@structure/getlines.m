function l=getlines(s,index)
%STRUCTURE/GETLINES getlines(s,index) returns s's lines with index
%  Each line is represented by two homogeneous 3D points.
%  See also getpluckerlines.

if nargin==1,
  l=s.lines;
else
  l=s.lines(:,index);
end



