function i=changecsystem(i,T);
% IMAGEDATA/changecsystem i=changecsystem(i,T) changes coordinate system
% T is a 3x3 projective transformation

if ~isempty(i.points)
  i.points = T * i.points;
end
if ~isempty(i.lines)
  if size(i.lines,1)==3,
    i.lines = inv(T)' * i.lines;
  else
    i.lines = blkdiag(T,T) * i.lines;
  end
end

for q=1:size(i.conics,2);
  m=v2m(i.conics(:,q));
  mnew = T*m*T';
  i.conics(:,q)=m2v(mnew);
end
