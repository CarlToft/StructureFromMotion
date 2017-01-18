function [pos2,residuals]=findpatchepi(im1,im2,pos1,F);
% blubb, 
patchsize=5;
stdev=1.5;
borderignore=12+patchsize;
[imsz1,imsz2]=size(im2);

a=stdev;
[x,y]=meshgrid(-3:3,-3:3);
filter=exp(- ((x).^2 + (y).^2)/a^2)/(a^2*pi);
im1filter=conv2(im1,filter,'same');
im2filter=conv2(im2,filter,'same');

nbrpoints=size(pos1,2);
pos2=NaN*ones(3,nbrpoints);
residuals=Inf*ones(1,nbrpoints);

for pp=1:nbrpoints,
  pt=pos1(:,pp);

%%%%%%%%%%%%%%%%

epiline=F'*pt;
epiline=epiline/norm(epiline(1:2)); 

horizontal=abs(epiline(1))<abs(epiline(2));

tt=[epiline(2);-epiline(1)];
if horizontal,
    xx=borderignore;
    yy=(-epiline(3)-epiline(1)*xx)/epiline(2);
    tt=tt*sign(tt(1));
else
    yy=borderignore;
    xx=(-epiline(3)-epiline(2)*yy)/epiline(1);
    tt=tt*sign(tt(2));
end

patch=im1filter(round(pt(2))+(-patchsize:patchsize),...
	  round(pt(1))+(-patchsize:patchsize));
sumpatch=sum(patch(:));

bestresidual=inf;
for nn=-5:1:5;
   uu=[xx;yy]+epiline(1:2)*nn;
   reslist=[];
   while (horizontal & uu(1)<imsz2-borderignore) | ...
        (~horizontal & uu(2)<imsz1-borderignore),

    if uu(1)>borderignore & uu(1)<imsz2-borderignore & ...
       uu(2)>borderignore & uu(2)<imsz1-borderignore,


      dpatch=im2filter(round(uu(2))+(-patchsize:patchsize),...
	         round(uu(1))+(-patchsize:patchsize));

      sumdpatch=sum(dpatch(:));if sumdpatch==0,sumdpatch=1;end
      tmp=sumpatch/sumdpatch*dpatch-patch;
      residual=sum(tmp(:).^2);
    else
      residual=Inf;
    end
    reslist=[reslist,residual];
    uu=uu+tt;
   end
   [residual,mindex]=min(reslist);
   if residual<bestresidual,
     destpos=[xx;yy]+(mindex-1)*tt+epiline(1:2)*nn;
     bestresidual=residual;
   end
end
if isfinite(bestresidual),
  [destpos,residual]=affineoptimize(im1,im2,pt,destpos,patchsize,stdev);

  residuals(pp)=residual;
  pos2(:,pp)=destpos;
end

%%%%%%%%%%%%%%%%%%%%


end %end of point-loop





