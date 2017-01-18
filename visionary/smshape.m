function [s,m]=smshape(imseq);
% [s,m]=smshape(imseq) solves the structure and motion
% problem for points with a shape based factorisation method.
% INPUT:
%   imseq - cell list containing imagedata objects
% OUTPUT:
%   s - structure
%   m - motion
% NB: All points must be visible in all images

nbrimages=length(imseq);
p0=getcommonfeatures(imseq);
nbrpoints=length(p0);

if nbrpoints<6,
  error('Too few points visible in all images');
end

udata = zeros(nbrpoints,3*nbrimages);
for i=1:nbrimages;
 u=getpoints(imseq{i},p0)';
 [u,s,v]=svd(u,0);
 udata(:,(3*i-2):3*i)=u;
end;

E=eye(nbrpoints);

for ii=1:30;
 %% Adjust U
 [us,ss,vs]=svd(udata,0);
 U=us(:,1:4);
 %% Adjust depth
 for jj=1:nbrimages;
  ui=udata(:,(3*jj-2):3*jj);
  Qu=E-U*U';
  M=[];
  for i=1:size(ui,2);
   M=[M;Qu*diag(ui(:,i))];
  end;
  [us,ss,vs]=svd(M'*M,0);
  q=vs(:,nbrpoints);
  q=q/sum(q);
  q=sqrt(size(q,1))*q/norm(q);

  ui=diag(q)*ui;
  [ui,ss,vs]=svd(ui,0);

  udata(:,(3*jj-2):3*jj)=ui;
 end;
end;

U=U';
%WARNING This has to be fixed properly!
U=pflat(U([2 3 4 1],:));
tmp=NaN*ones(4,size(imseq{1},1));
tmp(:,p0)=U;
s=structure(tmp);
m=reseclinear(imseq,s);
s=intseclinear(m,imseq,[],s); %intersect remaining points and other features

% end of main function



