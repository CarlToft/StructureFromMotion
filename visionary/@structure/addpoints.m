function s=addpoints(s,points)
%STRUCTURE/ADDPOINTS s=addpoints(s,points) adds homog. or cart. points

if size(points,1) == 3,
  s.points = [s.points,[points;ones(1,size(points,2))]];
else
  s.points = [s.points,points];
end

