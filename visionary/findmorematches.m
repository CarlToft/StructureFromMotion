function newimseq=findmorematches(images,imh,F12,F23,F13,threshold,intensitythreshold)
% blubb, 
patchsize=5;
stdev=1.5;
borderignore=12+patchsize;

nbrcorners=size(imh,1);
[imsz1,imsz2]=size(images{1});
newpts1=[];newpts2=[];newpts3=[];
for i=1:nbrcorners,
  pt=getpoints(imh,i);  
%figure(2);clf;plot(imagedata(images{2},pt));
%figure(1);clf;plot(imagedata(images{1}));

  epiline1=F12*pt;
  epiline1=epiline1/norm(epiline1(1:2)); 

  horizontal=abs(epiline1(1))<abs(epiline1(2));

  tt=[epiline1(2);-epiline1(1)];
  if horizontal,
    xx=borderignore;
    yy=(-epiline1(3)-epiline1(1)*xx)/epiline1(2);
    tt=tt*sign(tt(1));
  else
    yy=borderignore;
    xx=(-epiline1(3)-epiline1(2)*yy)/epiline1(1);
    tt=tt*sign(tt(2));
  end

  patch=double(images{2}(round(pt(2))+(-patchsize:patchsize),...
	      round(pt(1))+(-patchsize:patchsize)));
  sumpatch=sum(patch(:));

  bestresidual=inf;
  for nn=-8:1:8;
   uu=[xx;yy]+epiline1(1:2)*nn;
   reslist=[];
   while (horizontal & uu(1)<imsz2-borderignore) | ...
        (~horizontal & uu(2)<imsz1-borderignore),

    if uu(1)>borderignore & uu(1)<imsz2-borderignore & ...
       uu(2)>borderignore & uu(2)<imsz1-borderignore,


      dpatch=double(images{1}(round(uu(2))+(-patchsize:patchsize),...
	      round(uu(1))+(-patchsize:patchsize)));

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
     destpos1=[xx;yy]+(mindex-1)*tt+epiline1(1:2)*nn;
     bestresidual=residual;
   end
  end
  if ~isinf(bestresidual),
    [destpos1,residual]=affineoptimize(double(images{2}),double(images{1}),pt,destpos1,patchsize,stdev);
  else
    residual=Inf;
  end
  %figure(1);clf;plot(imagedata(images{1},destpos1,epiline1));zoom on;residual
  if residual<intensitythreshold, %image 3 search

    epiline3=F23'*pt;
    epiline3=epiline3/norm(epiline3(1:2));
    epiline3b=F13'*destpos1;
    epiline3b=epiline3b/norm(epiline3b(1:2));

    [u,s,v]=svd([epiline3,epiline3b]');
    destpos3=pflat(v(:,3));
    if destpos3(1)>borderignore & destpos3(1)<imsz2-borderignore & ...
           destpos3(2)>borderignore & destpos3(2)<imsz1-borderignore,
      [destpos3,residual]=affineoptimize(double(images{2}),double(images{3}),pt,destpos3,patchsize,stdev);
    else
      residual=inf;
    end
    if residual<intensitythreshold,
        newpts1=[newpts1,destpos1];
        newpts2=[newpts2,pt];
        newpts3=[newpts3,destpos3];
    else
%%%%% not found directly, look further

      horizontal=abs(epiline3(1))<abs(epiline3(2));

      tt=[epiline3(2);-epiline3(1)];
      if horizontal,
        xx=borderignore;
        yy=(-epiline3(3)-epiline3(1)*xx)/epiline3(2);
        tt=tt*sign(tt(1));
      else
        yy=borderignore;
        xx=(-epiline3(3)-epiline3(2)*yy)/epiline3(1);
        tt=tt*sign(tt(2));
      end

      bestresidual=inf;
      for nn=-8:1:8,
       uu=[xx;yy]+epiline3(1:2)*nn;reslist=[];
       while (horizontal & uu(1)<imsz2-borderignore) | ...
            (~horizontal & uu(2)<imsz1-borderignore),
        if uu(1)>borderignore & uu(1)<imsz2-borderignore & ...
           uu(2)>borderignore & uu(2)<imsz1-borderignore,
%          abs(epiline3b'*[uu;1])<5*threshold,
%          sqrt(sum((transpt(1:2)-uu(1:2)).^2))<4*threshold,

          dpatch=double(images{3}(round(uu(2))+(-patchsize:patchsize),...
  	      round(uu(1))+(-patchsize:patchsize)));

          sumdpatch=sum(dpatch(:));if sumdpatch==0,sumdpatch=1;end
          gscale=sumpatch/sumdpatch;
          tmp=gscale*dpatch-patch;
          residual=sum(tmp(:).^2);
        else
          residual=Inf;
        end

        reslist=[reslist,residual];
        uu=uu+tt;
       end %while
       [residual,mindex]=min(reslist);
       if residual<bestresidual,
         destpos3=[xx;yy]+(mindex-1)*tt+epiline3(1:2)*nn;
         bestresidual=residual;
       end
      end %for nn

      if ~isinf(bestresidual),
        [destpos3,residual]=affineoptimize(double(images{2}),double(images{3}),pt,destpos3,patchsize,stdev);
        if residual<intensitythreshold & abs(epiline3b'*destpos3)<8*threshold,
          newpts1=[newpts1,destpos1];
          newpts2=[newpts2,pt];
          newpts3=[newpts3,destpos3];
        end
      end %isnan(bestresidual)

    end %else, not found directly ->search
  end %image 3 search along epipolar
end

newimseq={imagedata(sparse(imsz1,imsz2),newpts1),imagedata(sparse(imsz1,imsz2),newpts2),imagedata(sparse(imsz1,imsz2),newpts3)};

%qq=1;figure(qq);clf;plot(imagedata(double(images{qq}),destpos1,epiline1));zoom on;hold on;
%qq=2;figure(qq);plot(imagedata(double(images{qq}),pt));zoom on;
%qq=3;figure(qq);plot(imagedata(double(images{qq}),destpos3,[epiline3,epiline3b]));zoom on;
