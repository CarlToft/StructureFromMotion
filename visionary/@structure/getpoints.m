function p=getpoints(s,index)
%STRUCTURE/GETPOINTS getpoints(s,index) returns s's points with index
%  in homog. coordinates.

if nargin==1,
  p=s.points;
else
  p=s.points(:,index);
end



