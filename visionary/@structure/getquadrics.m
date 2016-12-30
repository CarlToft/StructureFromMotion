function q=getquadrics(s,index)
%STRUCTURE/GETQUADRICS getquadrics(s,index) returns s's quadrics with index
%  in dual coordinates.

if nargin==1,
  q=s.quadrics;
else
  q=s.quadrics(:,index);
end



