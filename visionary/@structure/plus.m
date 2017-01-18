function s=plus(s1,s2)
%STRUCTURE/PLUS s=plus(s1,s2) adds s2's features to s1
s = structure(s1);
s2 = structure(s2);

s.points = [s.points,getpoints(s2)];
s.lines = [s.lines,getlines(s2)];
s.quadrics = [s.quadrics,getquadrics(s2)];

