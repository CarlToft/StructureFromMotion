function [str,mot,imseq,imseq0]=randomscene(nbrimages,nbrpoints,nbrlines,nbrconics,option)
% RANDOMSCENE creates a random scene and images with random camera positions
% [str,mot,imseq,imseq0]=randomscene(nbrimages,nbrpoints,nbrlines,nbrconics,option)
% INPUT:
%   nbrimages - number of images
%   nbrpoints - number of points
%   nbrlines - number of lines
%   nbrconics - numbder of conics
%   option:
%    'affine' - use the affine camera model
%    'constant' - constant intrinsic parameters
%    'zeroprincipal' - principal point at zero
%    'skew' - varying intrinsic parameters, except skew is zero
%    'distance=X' - set nominal camera distance to 'X'
%    'cube=X' - set nominal cube to [-X,X]
%    'noise=X' - add noise to image-points with std 'X'
%    'calibrated' - calibration matrices equal to identity
%    'norotation' - no relative rotation, hence translating cameras
%    'conics' - Quadrics are rank 3 disk-quadrics, i.e. conics in a plane
%    'closecameras' - cameras are close to previous camera
%    'planar' - points will be on z=0-plane
% OUTPUT:
%   str - structure
%   mot - motion
%   imseq - cell list of imagedata objects
%   imseq0 - unperturbated version of imseq
%    
% Default option is perspective camera with skew 0 and aspect ratio 1,
% with focal lengths around 1000, and 3D-features within the cube [-500,500]
% and camera distance around 1000 to origo.

if nargin<1,
  nbrimages=5;
end
if nargin<2,
  nbrpoints=10;
end
if nargin<3,
  nbrlines=0;
end
if nargin<4,
  nbrconics=0;
end

planar=0;
distance = 1000;
cube = 500;
noise = 0;
norotation=0;
affinecamera=0;
linesampling=10;
conicsampling=10;
redo=0;
camidentity=0;
diskquads=0;
zeroprincipal=0;
constantinternals=0;

closecameras=0;

if nargin>=5, %option set?
  if strmatch('planar',option);
    planar=1;
  end
  if strmatch('closecameras',option);
    closecameras=1;
  end
  if strmatch('affine',option);
    affinecamera=1;
  end
  if strmatch('constant',option);
    constantinternals=1;
  end
  if strmatch('norotation',option);
    norotation=1;
  end
  if strmatch('zeroprincipal',option);
    zeroprincipal=1;
  end

  if strmatch('skew',option);
    disp('skew option not implemented. PLEASE do!');
    keyboard;
  end
  if strmatch('distance',option);
    q=strmatch('distance=',option);
    strdist = option{q}(10:length(option{q}));
    distance=str2num(strdist);
  end
  if strmatch('cube',option);
    q=strmatch('cube=',option);
    strcube = option{q}(6:length(option{q}));
    cube=str2num(strcube);
  end
  if strmatch('noise',option);
    q=strmatch('noise=',option);
    strnoise = option{q}(7:length(option{q}));
    noise=str2num(strnoise);
  end
  if strmatch('calibrated',option);
    camidentity=1;
  end
  if strmatch('conics',option);
    diskquads=1;
  end
end

%create random points
point3d=[2*cube*(rand(3,nbrpoints)-0.5);ones(1,nbrpoints)];
if planar,
    point3d(3,:)=0;
end

%create random lines
cont=1;
while cont==1,
  line3d=[2*cube*(rand(3,nbrlines)-0.5);ones(1,nbrlines);...
	2*cube*(rand(3,nbrlines)-0.5);ones(1,nbrlines)];
  if nbrlines == 0,
    cont=0;
  else
    tmp=sqrt(min(sum((line3d(1:4,:)-line3d(5:8,:)).^2)));
    if tmp>cube/10,
      cont=0;
    end
  end
end

% create random 3D-ellipsoids
quadric = zeros(10,nbrconics);
for i=1:nbrconics,
  % axes
  Q=diag([(1/750+1/50*rand(3,1)).^2;-1]);

  % center
  pos=2*cube*(rand(3,1)-0.5);

  % rotation
  angles=rand(3,1)*2*pi;
  rot=getrotation(angles(1),angles(2),angles(3));
  T=[rot,pos;0 0 0 1];

  invT=inv(T);
  Q=invT'*Q*invT; % locus

  if diskquads==0,
    Ql=m2v(inv(Q));
  else
%   [u,d,v]=svd(inv(Q));
%   d(2,2)=0;
%   Ql=m2v(u*d*v');
   [u,t]=schur(Q);tmp=1./diag(t);tmp(1)=0;
   Ql=m2v(u*diag(tmp)*u');
  end
  quadric(:,i)=Ql;
end

mot=motion;
str=structure(point3d,line3d,quadric);
imseq={};
imseq0={};

% generate camera positions
for l=1:nbrimages
  t=distance*[randn(1,1)*0.1 randn(1,1)*0.1 (1+randn(1,1)*0.1)]';
  if norotation,
      R=eye(3);
  elseif l>1 & closecameras,
      R=R*expm(skew(randn(3,1)/5));
  else
      R=randrot;
  end
  if affinecamera,
    R(3,:)=0;
  end
  f=1000+randn*50;
  if camidentity,
    K=eye(3);
  elseif zeroprincipal
    K=[f 0 0;
     0 f 0;
     0 0 1];
  else
    K=[f 0 100*randn;
     0 f 100*randn;
     0 0 1];
  end

  if constantinternals,
    if l==1,
      Kfirst=K;
    else
      K=Kfirst;
    end
  end
  pcam=K*[R t];
  if affinecamera,
    pcam=pcam/pcam(3,4);
    [tmp1,tmp2]=rq(pcam(1:2,1:2));
    if det(tmp2)==-1,
      [tmp3,tmp4]=eig(tmp2);
      pcam(1:2,1:2)=tmp1*tmp3;
    end
  end
  mot=addcameras(mot,pcam);
  im = project(str,motion({pcam}));
  imseq0={imseq0{:},im};

  % perturbate points
  pts = pflat(getpoints(im));
  im = clearpoints(im);
  if nbrpoints>0,
    pts(1:2,:)=pts(1:2,:)+noise*randn(2,size(pts,2));
  end
  im = addpoints(im,pts);

  % perturbate lines
  lend1=pflat(pcam*line3d(1:4,:));
  lend2=pflat(pcam*line3d(5:8,:));
  nline=zeros(3,nbrlines);
  for j=1:nbrlines;
   d=lend2(:,j)-lend1(:,j);
   normal=normr([d(2), -d(1)]);

   lx=lend1(1,j):d(1)/(linesampling-1)-sign(d(1))*eps:lend2(1,j);
   ly=lend1(2,j):d(2)/(linesampling-1)-sign(d(2))*eps:lend2(2,j);

   linenoise=noise*randn(size(lx));
   if size(lx,2)==size(ly,2),
     lx=lx+normal(1)*linenoise;
     ly=ly+normal(2)*linenoise;
     if noise ==0,
       [u,s,v]=svd([lx',ly',ones(size(lx'))]);
       nline(:,j)=v(:,3);
     else
       nline(:,j)=gethomogeneouslines(fitline([lx;ly;ones(1,size(lx,2))]));
     end
   else
     redo = 1; %a line is projected as "point" ->failed
   end
  end
  im = clearlines(im);
  im = addlines(im,nline);

  % perturbate conics
  if nbrconics>0,
    co=getconics(im);

    for iq=1:nbrconics,
      c = inv(v2m(co(:,iq)));
      pts=getpointsonconic([c(1,1),c(1,2)*2,c(2,2),...
           c(1,3)*2,c(2,3)*2,c(3,3)],conicsampling);

tmp=v2m(co(:,iq));
tmp(1:2,3)=0;tmp(3,1:2)=0;tmp=inv(tmp);
     if isempty(pts) | min(svd(tmp(1:2,1:2)/tmp(3,3)))<5e-7,
        redo=1;
      else
        pts=pts+noise*randn(2,size(pts,2));
        tmp=fitconic([pts;ones(1,size(pts,2))]);
        if isempty(tmp),
	  redo=1;
	else
          co(:,iq)=pflat(getconics(tmp));

tmp=v2m(co(:,iq));
tmp(1:2,3)=0;tmp(3,1:2)=0;tmp=inv(tmp);
if isempty(pts) | min(svd(tmp(1:2,1:2)/tmp(3,3)))<5e-7,
redo=1;
end
	end
%        [u,s,v]=svd([pts(1,:).^2;2*pts(1,:).*pts(2,:);pts(2,:).^2;...
%                   2*pts(1,:);2*pts(2,:);ones(1,size(pts,2))]');
%        co(:,iq)=pflat(m2v(inv(v2m(v(:,6))))));
       end
    end
    im = clearconics(im);
    im = addconics(im,co);
  end

  imseq={imseq{:},im};
end;


if redo==1, %failed somehow
  [str,mot,imseq,imseq0]=randomscene(nbrimages,nbrpoints,nbrlines,...
      nbrconics,option);
end



