function [l,lc]=gethomogeneouslines(i,index)
%IMAGEDATA/GETHOMOGENEOUSLINES [l,lc]=gethomogeneouslines(i,index) returns i's lines l with index
%  in homog. coordinates and covariances lc

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

if size(l,1)==6,
  ltmp=l;
  l=zeros(3,size(ltmp,2));
  for ii=1:size(ltmp,2);
    l(:,ii)=cross(ltmp(1:3,ii),ltmp(4:6,ii));
  end
end


