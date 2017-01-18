function sz = size(i,index)
%IMAGEDATA/SIZE sz = size(i,index) returns nbr of pts, lines, conics

sz=[size(i.points,2),size(i.lines,2),size(i.conics,2)];

if nargin==2,
  sz = sz(index);
end
