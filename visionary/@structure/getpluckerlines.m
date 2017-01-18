function l=getpluckerlines(s,index)
%STRUCTURE/GETPLUCKERLINES getpluckerlines(s,index) returns s's
% lines (with index) in plucker coordinates.

if nargin==1,
  tmp=s.lines;
else
  tmp=s.lines(:,index);
end

l=zeros(6,size(tmp,2));
for i=1:size(tmp,2);
  tt=[tmp(1:4,i)';tmp(5:8,i)'];
  l(:,i)=[det(tt(:,[1,2]));det(tt(:,[1,3]));det(tt(:,[1,4]));...
	  det(tt(:,[3,4]));det(tt(:,[4,2]));det(tt(:,[2,3]))];
end

l=normc(l);





