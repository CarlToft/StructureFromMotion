function s=intsecpoints(m,imseq,points,option);
% INTSECPOINTS s=intsecpoints(m,imseq,points,option) calculates 3D points
% with linear method
% INPUT:
%   m - motion object
%   imseq - cell array of imagedata objects
%   points - (optional) specifies points to be intersected. Otherwise all points.
%            if 'points' is a structure object, already existing points are copied
%   option: 'nocoordtrans' - No coordinate transformation
% OUTPUT:
%   s - structure object containing reconstructed 3D points
% The results may be inaccurate.

if nargin<=2 | isempty(points)
  points = 1:size(imseq{1},1);
end

if nargin<=3
  option = [];
end
if isa(points,'structure');
  nbrpoints = size(imseq{1},1);
  Upts = NaN*ones(4,nbrpoints);
  if size(points,1)>0,
    Upts(:,1:size(points,1)) = getpoints(points);
  end
  points = 1:nbrpoints;
else
  nbrpoints=length(points);
  Upts = NaN*ones(4,nbrpoints);
end
Tnorm=cell(1,length(imseq));
for i=1:length(imseq),
  if isempty(strmatch('nocoordtrans',option)),
    Tnorm{i}=getnormtrans(imseq{i});
  else
    Tnorm{i}=eye(3);
  end
end
for qq=1:nbrpoints;
 if ~isfinite(Upts(1,qq)), %already existing?
  M=[];
  jj=0;
  for i=1:length(imseq);
   P=getcameras(m,i);
   u=getpoints(imseq{i},points(qq));
   if isfinite(P(1,1)) & isfinite(u(1)),
     %normalise
     P=Tnorm{i}*P;
     u=Tnorm{i}*u;
     jj=jj+1;
     M=[M zeros(3*(jj-1),1); ...
       P/norm(P,'fro') zeros(3,jj-1) psphere(u)];
   end;
  end;
  if jj>1, %too few views?
   [u,s,v]=svd(M);
   [U,alpha]=psphere(v(1:4,4+jj));
   dd=-v(5:(4+jj),4+jj)/alpha;
   if sum(sign(dd))<0;
     dd=-dd;
     U=-U;
   end;
   Upts(:,qq)=U;
  end
 end %already existing

end

s=structure(Upts);
