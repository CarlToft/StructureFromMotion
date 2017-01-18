function [sut,mut,step,lambda,oldf,newf,resnew,didbetter]=bundlestepplc(im,s,m,lambda,mode,caliboptions,lockcsystem,pp,focal,ar,zplane);


%calculate indices
[mi,si,imi,offsetmotion,offsetsm,offsetres]=calcindexplc(im,s,m,mode,zplane);

%calculate residual and gradient
[res,resgrad]=calcresgradplc(im,s,m,imi,si,mi,offsetres,offsetsm,pp,focal,ar,zplane);

oldf=norm(res);
resgrad=sparse(resgrad);
if mode == 1 | mode == 4 | mode ==5 | mode==7, %ordinary projective or only
                                     %structure or only motion or
                                     %knownrotation
				     
%sum(log(svd(full(resgrad)))<-10)-length(mi)-15
%keyboard;

  dx=-((resgrad'*resgrad+lambda*speye(offsetsm)))\(resgrad'*res);
  [sut,mut]=dolocparamplc(s,m,si,mi,dx,zplane);
elseif mode == 2, %affine
  nbrvar=size(resgrad,2);
  nbrcam=size(mi,1);
  nbrstr=nbrvar-12*nbrcam;
  dPdX=sparse(nbrvar,9*nbrcam+nbrstr);
  index=(9*nbrcam+1):(9*nbrcam+nbrstr);
  dPdX(3*nbrcam+index,index)=speye(nbrstr);

  dpdx=sparse(12,9);
  dpdx(1,1)=1;dpdx(2,2)=1;
  dpdx(4,3)=1;dpdx(5,4)=1;
  dpdx(7,5)=1;dpdx(8,6)=1;
  dpdx(10,7)=1;dpdx(11,8)=1;dpdx(12,9)=1;
  for qq=1:nbrcam;
    dPdX(12*qq-11:12*qq,9*qq-8:9*qq)=dpdx;
  end

  resgrad2=resgrad*dPdX;
  dy=-((resgrad2'*resgrad2+lambda*speye(9*nbrcam+nbrstr)))\(resgrad2'*res); 
  dx=zeros(offsetsm,1);
  dx(12*nbrcam+1:offsetsm)=dy(index);
  for qq=1:nbrcam;
    q=9*qq-8;
    dx(qq*12-11:qq*12)=...
      [dy(q:q+1);0;dy(q+2:q+3);0;dy(q+4:q+5);0;dy(q+6:q+8)];
  end
  [sut,mut]=dolocparamplc(s,m,si,mi,dx,zplane);
elseif mode == 3 | mode == 6, %autocalibration
  % Extract K,R and t from camera matrices
  [K,R,t]=ackrtextract(m);
  % Generate dPDX derivatives of all parameters with respect
  % to fewer (incorporating autocalibration constraints).
  nbrvar=size(resgrad,2);
  nbrcam=size(mi,1);
  nbrstr=nbrvar-12*nbrcam;
  if lockcsystem==0,
    nbrcamfreedom = 6*nbrcam; 
  else
    nbrcamfreedom = 6*(nbrcam-1)-1; %first camera is fix
  end
  nbrcamfreedom = nbrcamfreedom +sum(caliboptions==2);
  nbrcamfreedom = nbrcamfreedom +sum(caliboptions==3)*nbrcam;
  %
  dPdX=sparse(nbrvar,nbrcamfreedom+nbrstr);
  index1=(12*nbrcam+1):(12*nbrcam+nbrstr);
  index2=(nbrcamfreedom+1):(nbrcamfreedom+nbrstr);
  dPdX(index1,index2)=speye(nbrstr);
  %
  [dPdx,ppout,focalout,arout] = accalcdPdx(K,R,t,caliboptions,lockcsystem,pp,focal,ar); % <<!! Här räknas derivatorna
    % av alla kameramatriselement, m a p våra nya kameravariabler ut.
  %
  dPdX(1:(12*nbrcam),1:nbrcamfreedom)=dPdx;

  % Generate the new gradient matrix
  resgrad2=resgrad*dPdX;
  if ~isempty(focal),
    [focalm,focaln]=size(focalout);
    resgrad2(end+[1:focalm],1:focaln)=focalout;
  end
  if ~isempty(ar),
    [arm,arn]=size(arout);
    resgrad2(end+[1:arm],1:arn)=arout;
  end
  if ~isempty(pp),
    [ppm,ppn]=size(ppout);
    resgrad2(end+[1:ppm],1:ppn)=ppout;
  end

%for fully calibrated:
%sum(log(svd(full(resgrad2)))<-10)-7
%keyboard;
  dy=-((resgrad2'*resgrad2+lambda*speye(size(resgrad2,2))))\(resgrad2'*res);

  [mut] = aclocparam(K,R,t,dy(1:nbrcamfreedom),caliboptions,lockcsystem);

  [sut,mskrap]=dolocparamplc(s,m,si,mi,[zeros(12*nbrcam,1);dy(index2)],zplane);
else
  error('Unknown mode');
end


[resnew]=calcresplc(im,sut,mut,imi,si,mi,offsetres,offsetsm,pp,focal,ar);

newf=norm(resnew);

if newf>oldf,
 ii=0;
 lambda=lambda*2;
 while ( (ii<10) & (newf>oldf) ),
  if mode ~= 3 & mode ~=6,
    dx=dx/2;
    [sut,mut]=dolocparamplc(s,m,si,mi,dx,zplane);
  else
    dy=dy/2;
    [mut] = aclocparam(K,R,t,dy(1:nbrcamfreedom),caliboptions,lockcsystem);
    [sut,mskrap]=dolocparamplc(s,m,si,mi,[zeros(12*nbrcam,1);dy(index2)],zplane);
  end;
  [resnew]=calcresplc(im,sut,mut,imi,si,mi,offsetres,offsetsm,pp,focal,ar);
  newf=norm(resnew);
  ii=ii+1;
 end;
else
  lambda=lambda/3;   
end;

if oldf<newf,
  sut=s;
  mut=m;
%  disp('ingen förbättring');
  resnew=res;
  didbetter=0;
else
  didbetter=1;
end;

if mode == 3 | mode ==6,
  step=norm(dy);
else
  step=norm(dx);
end

