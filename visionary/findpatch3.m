function [pos2,reslist]=findpatch3(im1,im2,pos1,pred2,searchsize)
% [pos2,residual]=findpatch3(im1,im2,pos1)
% Finds a patch in image im2,
% using a patch around pos1 in image im1
% by means of affine subpixel correlation
% Input:
%  im1,im2 - image matrices
%  pos1 - positions in image im1
%  pred2 - predicted position (optional)
%  searchsize - size of search area (optional)
% Output:
%  pos2 - positions in image im2
%  residual - intensity residual between patches
% See also: findpatch,findpatch2

if nargin<=3 | isempty(pred2),
  pred2=pos1;
end;
if nargin<=4,
  searchsize=100;
end;

patchsize=5;
stdev=1.5;
n=2*patchsize+1;

alfa1=2/180*pi;
alfa2=4/180*pi;
scale1=0.06;
scale2=0.12;

nbrpoints=size(pos1,2);
pos2=NaN*ones(3,nbrpoints);
reslist=Inf*ones(1,nbrpoints);

borderignore=8;

a=stdev;
[x,y]=meshgrid(-3:3,-3:3);
filter=exp(- ((x).^2 + (y).^2)/a^2)/(a^2*pi);
im1filter=conv2(im1,filter,'same');
im2filter=conv2(im2,filter,'same');

im2ones=conv2(im2filter,ones(n),'valid');
im2ones(find(im2ones==0))=1;
im2rho=1./im2ones;
im2squared=conv2(im2filter.^2,ones(n),'valid');



for pp=1:nbrpoints,

  pt1=pos1(:,pp);
  predpos=pred2(:,pp);


  if isnan(pt1(1)) | sum(pt1(1:2)<2*patchsize+1)>0 | ...
	pt1(1)>size(im1,2)-2*patchsize | ...
	pt1(2)>size(im1,1)-2*patchsize,
   residual=Inf;
   pt2=NaN*ones(3,1);
else %inside image
%%%%%%%%%%%%%%%%%%%




index1=max(1+borderignore,round(predpos(2))-searchsize):min(size(im2,1)-borderignore,round(predpos(2))+searchsize);
index2=max(1+borderignore,round(predpos(1))-searchsize):min(size(im2,2)-borderignore,round(predpos(1))+searchsize);
searchpatch=im2filter(index1,index2);

searchsquared=im2squared(index1(1+patchsize:end-patchsize)-patchsize,...
		 index2(1+patchsize:end-patchsize)-patchsize);

rho=im2rho(index1(1+patchsize:end-patchsize)-patchsize,...
		 index2(1+patchsize:end-patchsize)-patchsize);


rhosquared=rho.^2;

bestdist=Inf;
for method=1:5,

  switch method,
    case 1,
      patch=im1filter(round(pt1(2))+(-patchsize:patchsize),...
	      round(pt1(1))+(-patchsize:patchsize));
      theR=eye(2);
    case 2,
      tmp=imrotate(im1filter(round(pt1(2))+(-2*patchsize:2*patchsize),...
	round(pt1(1))+(-2*patchsize:2*patchsize)),alfa1*180/pi,'bicubic','crop');
      patch=tmp(patchsize+1:3*patchsize+1,patchsize+1:3*patchsize+1);
      theR=[cos(alfa1),sin(alfa1);-sin(alfa1),cos(alfa1)];
    case 3,
      tmp=imrotate(im1filter(round(pt1(2))+(-2*patchsize:2*patchsize),...
	  round(pt1(1))+(-2*patchsize:2*patchsize)),-alfa1*180/pi,'bicubic','crop');
      patch=tmp(patchsize+1:3*patchsize+1,patchsize+1:3*patchsize+1);
      theR=[cos(alfa1),sin(-alfa1);sin(alfa1),cos(alfa1)];
    case 4,
      tmp=imresize(im1filter(round(pt1(2))+(-patchsize:patchsize),...
	  round(pt1(1))+(-patchsize:patchsize)),1+scale1,'bicubic');
      iii=(size(tmp,1)+1)/2+(-patchsize:patchsize);
      patch=tmp(iii,iii);
      theR=diag(1+[scale1,scale1]);
    case 5,
      tmp=imresize(im1filter(round(pt1(2))+(-2-patchsize:patchsize+2),...
	  round(pt1(1))+(-2-patchsize:2+patchsize)),1-scale1,'bicubic');
      iii=(size(tmp,1))/2+(-patchsize:patchsize);
      patch=tmp(iii,iii);
      theR=diag(1-[scale1,scale1]);
    case 6,
      tmp=imrotate(im1filter(round(pt1(2))+(-2*patchsize:2*patchsize),...
	round(pt1(1))+(-2*patchsize:2*patchsize)),alfa2*180/pi,'bicubic','crop');
      patch=tmp(patchsize+1:3*patchsize+1,patchsize+1:3*patchsize+1);
      theR=[cos(alfa2),sin(alfa2);-sin(alfa2),cos(alfa2)];
    case 7,
      tmp=imrotate(im1filter(round(pt1(2))+(-2*patchsize:2*patchsize),...
	  round(pt1(1))+(-2*patchsize:2*patchsize)),-alfa2*180/pi,'bicubic','crop');
      patch=tmp(patchsize+1:3*patchsize+1,patchsize+1:3*patchsize+1);
      theR=[cos(alfa2),sin(-alfa2);sin(alfa2),cos(alfa2)];
    case 8,
      tmp=imresize(im1filter(round(pt1(2))+(-patchsize:patchsize),...
	  round(pt1(1))+(-patchsize:patchsize)),1+scale2,'bicubic');
      iii=size(tmp,1)/2+(-patchsize:patchsize);
      patch=tmp(iii,iii);
      theR=diag(1+[scale2,scale2]);
    case 9,
      tmp=imresize(im1filter(round(pt1(2))+(-3-patchsize:patchsize+3),...
	  round(pt1(1))+(-3-patchsize:3+patchsize)),1-scale2,'bicubic');
      iii=size(tmp,1)/2+(-patchsize:patchsize);
      patch=tmp(iii,iii);
      theR=diag(1-[scale2,scale2]);
  end

  patch = flipud(fliplr(patch));
  patchsquared=sum(patch(:).^2);if patchsquared==0,patchsquared=1;end;
  patchsum=sum(patch(:));

  search_vs_patch=conv2(searchpatch,patch,'valid');

  distmap=patchsum^2*(rho.^2).*searchsquared+patchsquared-...
            2*patchsum*rho.*search_vs_patch;
  dist=min(distmap(:));

  if dist<bestdist,
    [posxx,posyy]=find(distmap==dist);
    posxx=posxx(1)-1+index1(1)+patchsize;
    posyy=posyy(1)-1+index2(1)+patchsize;
%bestmethod=method;
    bestR=theR;
    bestdist=dist;
  end
end %loop over method

%bestmethod
destpos=[posyy;posxx];

[pt2,residual]=affineoptimize(im1,im2,pt1,destpos,patchsize,stdev,inv(theR));


%%%%%%%%%%%%%%%%%%%
end %insidee image


  reslist(pp)=residual;
  pos2(:,pp)=pt2;

end %end of point-loop

