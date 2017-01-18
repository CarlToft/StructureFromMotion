function [c,cc]=getconics(i,index)
%IMAGEDATA/GETCONICS [c,cc]=getconics(i,index) returns i's conics with index
%  in dual coordinates and covariances cc

if nargin==1,
  c=i.conics;
  cc=i.coniccov;
else
  c=i.conics(:,index);
  if nargout==2,
    cc=i.coniccov(index);
    if size(cc,2)==1;
      cc=cc{1};
    end
  end
end



