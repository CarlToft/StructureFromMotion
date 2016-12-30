function [l,lc]=getlines(i,index)
%IMAGEDATA/GETLINES [l,lc]=getlines(i,index) returns i's lines l with index
%  in stored format and covariances lc

if nargin==1,
  l=i.lines;
  lc=i.linecov;
else
  l=i.lines(:,index);
  if nargout==2,
    if isempty(i.linecov),
      lc=[1 0 0;0 1 0];
    else
      lc=i.linecov(index);
      if size(lc,2)==1;
      lc=lc{1};
      end
    end
  end
end



