function sz = size(s,index)
%STRUCTURE/SIZE sz = size(s,index) returns nbr of pts, lines, quadrics

sz=[size(s.points,2),size(s.lines,2),size(s.quadrics,2)];

if nargin==2,
  sz = sz(index);
end
