function s=intseclines(m,imseq,lines);
% INTSECLINES s=intseclines(m,imseq,lines,option) calculates 3D lines
% with linear method
% INPUT:
%   m - motion object
%   imseq - cell array of imagedata objects
%   lines - (optional) specifies lines to be intersected. Otherwise all lines.
%            if 'lines' is a structure object, already existing lines are copied
% OUTPUT:
%   s - structure object containing reconstructed 3D lines
% The results may be inaccurate.

if nargin<=2 | isempty(lines)
  lines = 1:size(imseq{1},2);
end

if isa(lines,'structure');
  nbrlines = size(lines,2);
  if nbrlines>0,
    Ldata = gethomogeneouslines(lines);
  else
    nbrlines=size(imseq{1},2);
    Ldata = NaN*ones(8,nbrlines);
  end
else
  nbrlines=length(lines);
  Ldata = NaN*ones(8,nbrlines);
end

nbrimages=length(imseq);

for qq=1:nbrlines;
 if ~finite(Ldata(1,qq)), %already existing
  jj=0;
  M=[];
  for i=1:nbrimages;
   P=getcameras(m,i);
   l=gethomogeneouslines(imseq{i},qq);
   if finite(P(1)) & finite(l(1)),
     plane = P'*l;
     M = [M;psphere(plane)'];
     jj=jj+1;
   end
  end
  if jj>2, %too few views?
   [u,ss,v]=svd(M);
   Ldata(:,qq)=[v(:,3);v(:,4)];
  end
 end %already existing
end

s=structure([],Ldata);




