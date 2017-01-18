function s=changecsystem(s,T);
% STRUCTURE/changecsystem s=changecsystem(s,T) changes coordinate system
% T is a 4x4 projective transformation

if ~isempty(s.points)
  s.points = T * s.points;
end
if ~isempty(s.lines)
  s.lines = [T,zeros(4);zeros(4),T] * s.lines;
end

for i=1:size(s.quadrics,2);
  m=v2m(s.quadrics(:,i));
  mnew = T*m*T';
  s.quadrics(:,i)=m2v(mnew);
end
