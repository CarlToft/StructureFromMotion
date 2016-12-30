function [s,m]=smaffine(imseq,option);
% [s,m]=smaffine(imseq,option) solves the structure and motion with affine cameras
% for points, lines and conics with factorisation.
% NB: Only features that are visible in all images are used.
% INPUT:
%   imseq - cell list containing imagedata objects
%   option:
%     'conic' - the quadrics are disk-quadrics, i.e. conics
% OUTPUT:
%   s - structure
%   m - motion
%
% See also: smaffineclos

if nargin<=1,
  option ='';
end

nbrimages=length(imseq);
[p0,l0,c0]=getcommonfeatures(imseq);

nbrpoints=length(p0);
nbrlines=length(l0);
nbrconics=length(c0);
nbrpc=nbrpoints+nbrconics;

% Check for minimality conditions (NB: This is not done strictly)
if nbrpc<5,
 if nbrpc<3 | nbrlines < 3;
   error(['Too few visible features: ', ...
           num2str(nbrpoints),' points ', ...
	   num2str(nbrlines),' lines ', ...
           num2str(nbrconics),' conics']);
 end
end

trans=cell(1,nbrimages);
for i=1:nbrimages;
  if isempty(strmatch('nocoordtrans',option)),
    trans{i}=getnormtrans(imseq{i},'affine');
  else
    trans{i}=eye(3);
  end
end

Utilde=[];  %notation according to Tomasi-Kanade
Vtilde=[];
centroid = zeros(2,nbrimages);
%scale = ones(2,nbrimages);

for i=1:nbrimages,
  im=changecsystem(imseq{i},trans{i});

  %points
  if nbrpoints>0
    pts=pflat(getpoints(im,p0));
    pts = pts(1:2,:);
  end

  %conics
  if nbrconics>0,
    u = getconics(im,c0);
    u = pflat(u(4:6,:));
    if nbrpoints>0,
      pts = [pts,u(1:2,:)];
    else
      pts = u(1:2,:);
    end
  end

  %relative coordinates
  centroid(:,i)=mean(pts')';
  utilde = pts(1,:) - centroid(1,i);
  vtilde = pts(2,:) - centroid(2,i);
%  scale(:,i) = [max(abs(utilde));max(abs(vtilde))];
%  utilde = utilde/scale(1,i);
%  vtilde = vtilde/scale(2,i);

  %lines
  if nbrlines>0,
    uline = pflat(gethomogeneouslines(im,l0));
%    uline = diag([scale(:,i)',1])*pflat(gethomogeneouslines(im,l0));
    udir = normc([uline(2,:);-uline(1,:)]);
    utilde = [utilde,udir(1,:)];
    vtilde = [vtilde,udir(2,:)];
  end
  Utilde = [Utilde;utilde];
  Vtilde = [Vtilde;vtilde];
end
if nbrlines>0,
  %retrieve direction scale-factor
  depth=rescalelines(Utilde,Vtilde,nbrimages,nbrpc,nbrlines);
  index=nbrpc+[1:nbrlines];
  Utilde(:,index)=Utilde(:,index).*depth;
  Vtilde(:,index)=Vtilde(:,index).*depth;
end

W= [Utilde;Vtilde];
[O1,sv,O2]=svd(W);
sv3=sv(1:3,1:3);
Rtilde = O1(:,1:3)*sqrt(sv3);
Stilde = sqrt(sv3)*O2(:,1:3)';

m = motion;
for i=1:nbrimages,

%  T = diag(scale(:,i));
%  pcam = [ T*[Rtilde(i,:);Rtilde(i+nbrimages,:)],...
%           centroid(:,i);zeros(1,3),1 ];
%  pcam = inv(trans{i})*[ T*[Rtilde(i,:);Rtilde(i+nbrimages,:)],...
%           centroid(:,i);zeros(1,3),1 ];
  pcam = inv(trans{i})*[ [Rtilde(i,:);Rtilde(i+nbrimages,:)],...
           centroid(:,i);zeros(1,3),1 ];

  m = addcameras(m,pcam);
end
tmp=NaN*ones(4,size(imseq{1},1));
tmp(:,p0)=[Stilde(:,1:nbrpoints);ones(1,nbrpoints)];
s=structure(tmp);
s=intseclinear(imseq,m,option,s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of main-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function depth=rescalelines(Utilde,Vtilde,nbrimages,nbrpc,nbrlines)

M=zeros(15*nbrimages,nbrimages,nbrlines);
total=nbrpc+nbrlines;

for bildnr=1:nbrimages,
 ind1=bildnr;
 ind2=mod(ind1,nbrimages)+1;
 ind3=mod(ind2,nbrimages)+1;

 Utemp=Utilde([ind1,ind2,ind3],:);
 Vtemp=Vtilde([ind1,ind2,ind3],:);
 [A0,A1,A2]=tritensor20(Utemp,Vtemp,nbrpc,nbrlines);
 a=A0(1,1);b=A0(1,2);c=A0(1,3);
 d=A0(2,1);e=A0(2,2);f=A0(2,3);
 g=A1(1,1);h=A1(1,2);i=A1(1,3);
 j=A1(2,1);k=A1(2,2);l=A1(2,3);
 m=A2(1,1);n=A2(1,2);o=A2(1,3);
 p=A2(2,1);q=A2(2,2);r=A2(2,3);

 for linenr=1:nbrlines,
  u1=Utemp(1,nbrpc+linenr);
  u2=Vtemp(1,nbrpc+linenr);
  v1=Utemp(2,nbrpc+linenr);
  v2=Vtemp(2,nbrpc+linenr);
  w1=Utemp(3,nbrpc+linenr);
  w2=Vtemp(3,nbrpc+linenr);

  row=bildnr*15-14;
  col1=bildnr;
  col2=mod(col1,nbrimages)+1;
  col3=mod(col2,nbrimages)+1;

  M(row,col1,linenr)=a*h*u2*l-a*k*u2*i-d*h*u1*l+d*k*u1*i-g*b*u2*l+g*e*u1*l+g*k*c*u2-g*k* ...
u1*f+j*b*u2*i-j*e*u1*i-j*h*c*u2+j*h*u1*f;
  M(row,col2,linenr)=a*e*i*v2-a*e*v1*l-a*h*f*v2+a*k*f*v1-d*b*i*v2+d*b*v1*l+d*h*c*v2-d*k*c ...
*v1+g*b*f*v2-g*e*c*v2-j*b*f*v1+j*e*c*v1;

  M(row+1,col1,linenr)=a*h*u2*o-a*n*u2*i-d*h*u1*o+d*n*u1*i-g*b*u2*o+g*e*u1*o+g*n*c*u2-g*n* ...
u1*f+m*b*u2*i-m*e*u1*i-m*h*c*u2+m*h*u1*f;
  M(row+1,col2,linenr)=-a*e*v1*o+a*n*f*v1+d*b*v1*o-d*n*c*v1-m*b*f*v1+m*e*c*v1;
  M(row+1,col3,linenr)=a*e*i*w1-a*h*f*w1-d*b*i*w1+d*h*c*w1+g*b*f*w1-g*e*c*w1;

  M(row+2,col1,linenr)=a*h*u2*r-a*q*u2*i-d*h*u1*r+d*q*u1*i-g*b*u2*r+g*e*u1*r+g*q*c*u2-g*q* ...
u1*f+p*b*u2*i-p*e*u1*i-p*h*c*u2+p*h*u1*f;
  M(row+2,col2,linenr)=-a*e*v1*r+a*q*f*v1+d*b*v1*r-d*q*c*v1-p*b*f*v1+p*e*c*v1;
  M(row+2,col3,linenr)=a*e*i*w2-a*h*f*w2-d*b*i*w2+d*h*c*w2+g*b*f*w2-g*e*c*w2;

  M(row+3,col1,linenr)=a*k*u2*o-a*n*u2*l-d*k*u1*o+d*n*u1*l-j*b*u2*o+j*e*u1*o+j*n*c*u2-j*n* ...
u1*f+m*b*u2*l-m*e*u1*l-m*k*c*u2+m*k*u1*f;
  M(row+3,col2,linenr)=-a*e*v2*o+a*n*f*v2+d*b*v2*o-d*n*c*v2-m*b*f*v2+m*e*c*v2;
  M(row+3,col3,linenr)=a*e*l*w1-a*k*f*w1-d*b*l*w1+d*k*c*w1+j*b*f*w1-j*e*c*w1;

  M(row+4,col1,linenr)=a*k*u2*r-a*q*u2*l-d*k*u1*r+d*q*u1*l-j*b*u2*r+j*e*u1*r+j*q*c*u2-j*q* ...
u1*f+p*b*u2*l-p*e*u1*l-p*k*c*u2+p*k*u1*f;
  M(row+4,col2,linenr)=-a*e*v2*r+a*q*f*v2+d*b*v2*r-d*q*c*v2-p*b*f*v2+p*e*c*v2;
  M(row+4,col3,linenr)=a*e*l*w2-a*k*f*w2-d*b*l*w2+d*k*c*w2+j*b*f*w2-j*e*c*w2;

  M(row+5,col1,linenr)=a*n*u2*r-a*q*u2*o-d*n*u1*r+d*q*u1*o-m*b*u2*r+m*e*u1*r+m*q*c*u2-m*q* ...
u1*f+p*b*u2*o-p*e*u1*o-p*n*c*u2+p*n*u1*f;
  M(row+5,col3,linenr)=a*e*o*w2-a*e*w1*r-a*n*f*w2+a*q*f*w1-d*b*o*w2+d*b*w1*r+d*n*c*w2-d*q*c ...
*w1+m*b*f*w2-m*e*c*w2-p*b*f*w1+p*e*c*w1;

  M(row+6,col1,linenr)=-g*k*u1*o+g*n*u1*l+j*h*u1*o-j*n*u1*i-m*h*u1*l+m*k*u1*i;
  M(row+6,col2,linenr)=-a*h*v2*o+a*k*v1*o+a*n*i*v2-a*n*v1*l+g*b*v2*o-g*n*c*v2-j*b*v1*o+j*n* ...
c*v1-m*b*i*v2+m*b*v1*l+m*h*c*v2-m*k*c*v1;
  M(row+6,col3,linenr)=a*h*l*w1-a*k*i*w1-g*b*l*w1+g*k*c*w1+j*b*i*w1-j*h*c*w1;

  M(row+7,col1,linenr)=-g*k*u1*r+g*q*u1*l+j*h*u1*r-j*q*u1*i-p*h*u1*l+p*k*u1*i;
  M(row+7,col2,linenr)=-a*h*v2*r+a*k*v1*r+a*q*i*v2-a*q*v1*l+g*b*v2*r-g*q*c*v2-j*b*v1*r+j*q* ...
c*v1-p*b*i*v2+p*b*v1*l+p*h*c*v2-p*k*c*v1;
  M(row+7,col3,linenr)=a*h*l*w2-a*k*i*w2-g*b*l*w2+g*k*c*w2+j*b*i*w2-j*h*c*w2;

  M(row+8,col1,linenr)=-g*n*u1*r+g*q*u1*o+m*h*u1*r-m*q*u1*i-p*h*u1*o+p*n*u1*i;
  M(row+8,col2,linenr)=a*n*v1*r-a*q*v1*o-m*b*v1*r+m*q*c*v1+p*b*v1*o-p*n*c*v1;
  M(row+8,col3,linenr)=a*h*o*w2-a*h*w1*r-a*n*i*w2+a*q*i*w1-g*b*o*w2+g*b*w1*r+g*n*c*w2-g*q*c ...
*w1+m*b*i*w2-m*h*c*w2-p*b*i*w1+p*h*c*w1;

  M(row+9,col1,linenr)=-j*n*u1*r+j*q*u1*o+m*k*u1*r-m*q*u1*l-p*k*u1*o+p*n*u1*l;
  M(row+9,col2,linenr)=a*n*v2*r-a*q*v2*o-m*b*v2*r+m*q*c*v2+p*b*v2*o-p*n*c*v2;
  M(row+9,col3,linenr)=a*k*o*w2-a*k*w1*r-a*n*l*w2+a*q*l*w1-j*b*o*w2+j*b*w1*r+j*n*c*w2-j*q*c ...
*w1+m*b*l*w2-m*k*c*w2-p*b*l*w1+p*k*c*w1;

  M(row+10,col1,linenr)=-g*k*u2*o+g*n*u2*l+j*h*u2*o-j*n*u2*i-m*h*u2*l+m*k*u2*i;
  M(row+10,col2,linenr)=-d*h*v2*o+d*k*v1*o+d*n*i*v2-d*n*v1*l+g*e*v2*o-g*n*f*v2-j*e*v1*o+j*n* ...
f*v1-m*e*i*v2+m*e*v1*l+m*h*f*v2-m*k*f*v1;
  M(row+10,col3,linenr)=d*h*l*w1-d*k*i*w1-g*e*l*w1+g*k*f*w1+j*e*i*w1-j*h*f*w1;

  M(row+11,col1,linenr)=-g*k*u2*r+g*q*u2*l+j*h*u2*r-j*q*u2*i-p*h*u2*l+p*k*u2*i;
  M(row+11,col2,linenr)=-d*h*v2*r+d*k*v1*r+d*q*i*v2-d*q*v1*l+g*e*v2*r-g*q*f*v2-j*e*v1*r+j*q* ...
f*v1-p*e*i*v2+p*e*v1*l+p*h*f*v2-p*k*f*v1;
  M(row+11,col3,linenr)=d*h*l*w2-d*k*i*w2-g*e*l*w2+g*k*f*w2+j*e*i*w2-j*h*f*w2;

  M(row+12,col1,linenr)=-g*n*u2*r+g*q*u2*o+m*h*u2*r-m*q*u2*i-p*h*u2*o+p*n*u2*i;
  M(row+12,col2,linenr)=d*n*v1*r-d*q*v1*o-m*e*v1*r+m*q*f*v1+p*e*v1*o-p*n*f*v1;
  M(row+12,col3,linenr)=d*h*o*w2-d*h*w1*r-d*n*i*w2+d*q*i*w1-g*e*o*w2+g*e*w1*r+g*n*f*w2-g*q*f ...
*w1+m*e*i*w2-m*h*f*w2-p*e*i*w1+p*h*f*w1;

  M(row+13,col1,linenr)=-j*n*u2*r+j*q*u2*o+m*k*u2*r-m*q*u2*l-p*k*u2*o+p*n*u2*l;
  M(row+13,col2,linenr)=d*n*v2*r-d*q*v2*o-m*e*v2*r+m*q*f*v2+p*e*v2*o-p*n*f*v2;
  M(row+13,col3,linenr)=d*k*o*w2-d*k*w1*r-d*n*l*w2+d*q*l*w1-j*e*o*w2+j*e*w1*r+j*n*f*w2-j*q*f ...
*w1+m*e*l*w2-m*k*f*w2-p*e*l*w1+p*k*f*w1;

  M(row+14,col2,linenr)=g*n*v2*r-g*q*v2*o-j*n*v1*r+j*q*v1*o-m*h*v2*r+m*k*v1*r+m*q*i*v2-m*q* ...
v1*l+p*h*v2*o-p*k*v1*o-p*n*i*v2+p*n*v1*l;
  M(row+14,col3,linenr)=g*k*o*w2-g*k*w1*r-g*n*l*w2+g*q*l*w1-j*h*o*w2+j*h*w1*r+j*n*i*w2-j*q*i ...
*w1+m*h*l*w2-m*k*i*w2-p*h*l*w1+p*k*i*w1;
 end
end

depth=zeros(nbrimages,nbrlines);
for linenr=1:nbrlines,
 [U,S,V]=svd(M(:,:,linenr));
 depth(:,linenr)=V(:,nbrimages);
end

function [A0,A1,A2]=tritensor20(Utilde,Vtilde,nbrpc,nbrlines);
%reconstuct A0,A1,A2 from three views

M=zeros(15*nbrpc+nbrlines,20);
R=zeros(15,20);

for i=1:nbrpc,
  x1=Utilde(1,i);
  x2=Vtilde(1,i);
  y1=Utilde(2,i);
  y2=Vtilde(2,i);
  z1=Utilde(3,i);
  z2=Vtilde(3,i);

      R(1,9) = y2;
      R(1,10) = -y1;
      R(1,13) = x2;
      R(1,14) = -x1;
      R(2,1) = x2;
      R(2,5) = -x1;
      R(2,9) = z1;
      R(2,11) = -y1;
      R(3,2) = x2;
      R(3,6) = -x1;
      R(3,9) = z2;
      R(3,12) = -y1;
      R(4,3) = x2;
      R(4,7) = -x1;
      R(4,10) = z1;
      R(4,11) = -y2;
      R(5,4) = x2;
      R(5,8) = -x1;
      R(5,10) = z2;
      R(5,12) = -y2;
      R(6,11) = z2;
      R(6,12) = -z1;
      R(6,17) = x2;
      R(6,18) = -x1;
      R(7,1) = -y2;
      R(7,3) = y1;
      R(7,13) = z1;
      R(7,15) = -x1;
      R(8,2) = -y2;
      R(8,4) = y1;
      R(8,13) = z2;
      R(8,16) = -x1;
      R(9,1) = z2;
      R(9,2) = -z1;
      R(9,17) = y1;
      R(9,19) = -x1;
      R(10,3) = z2;
      R(10,4) = -z1;
      R(10,17) = y2;
      R(10,20) = -x1;
      R(11,5) = -y2;
      R(11,7) = y1;
      R(11,14) = z1;
      R(11,15) = -x2;
      R(12,6) = -y2;
      R(12,8) = y1;
      R(12,14) = z2;
      R(12,16) = -x2;
      R(13,5) = z2;
      R(13,6) = -z1;
      R(13,18) = y1;
      R(13,19) = -x2;
      R(14,7) = z2;
      R(14,8) = -z1;
      R(14,18) = y2;
      R(14,20) = -x2;
      R(15,15) = z2;
      R(15,16) = -z1;
      R(15,19) = y2;
      R(15,20) = -y1;
  M(i*15-14:i*15,:)=R;
end

linec=zeros(1,8);
for i=1:nbrlines,
  u1=Utilde(1,nbrpc+i);
  u2=Vtilde(1,nbrpc+i);
  v1=Utilde(2,nbrpc+i);
  v2=Vtilde(2,nbrpc+i);
  w1=Utilde(3,nbrpc+i);
  w2=Vtilde(3,nbrpc+i);

      linec(1) = -u2*v2*w2;
      linec(2) = u2*v2*w1;
      linec(3) = u2*v1*w2;
      linec(4) = -u2*v1*w1;
      linec(5) = u1*v2*w2;
      linec(6) = -u1*v2*w1;
      linec(7) = -u1*v1*w2;
      linec(8) = u1*v1*w1;

  M(nbrpc*15+i,1:8)=linec;
end
[u,d,v]=svd(M);

param20=v(:,20)/v(1,20);

A0=zeros(2,3);
A1=zeros(2,3);
A2=zeros(2,3);
A0(1,1)=1;
A1(1,2)=1;
A2(1,3)=1;

A0(2,1)=param20(5);
A0(2,2)=param20(11);
A0(2,3)=-param20(9);
A1(2,1)=-param20(15);
A1(2,2)=param20(3);
A1(2,3)=param20(13);
A2(2,1)=param20(19);
A2(2,2)=-param20(17);
A2(2,3)=param20(2);

%change to good coordinate system
[u,s,v]=svd([A0;A1;A2]);
t=v*diag(1./diag(s));
A0=A0*t;
A1=A1*t;
A2=A2*t;


