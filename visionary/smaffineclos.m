function [s,m]=smaffineclos(imseq,option);
% [s,m]=smaffineclos(imseq,option) solves the structure and motion with affine cameras
% for points, lines and conics with closure constraints
% Missing data is handled automaticly
% INPUT:
%   imseq - cell list containing imagedata objects
%   option:
%     'conic' - the quadrics are disk-quadrics, i.e. conics
%     'nocoordtrans' - No coordinate transformation
% OUTPUT:
%   s - structure
%   m - motion
%
% See also: smaffine
if nargin<=1,
  option ='';
end

nbrimages=length(imseq);

if nbrimages<3,
  error('Too few images');
end

trans=cell(1,nbrimages);
for i=1:nbrimages;
  if isempty(strmatch('nocoordtrans',option)),
    trans{i}=getnormtrans(imseq{i},'affine');
  else
    trans{i}=eye(3);
  end
end
centers=zeros(6,0);
samecenter=zeros(1,nbrimages);
T=zeros(0,20);
i=0;
while i<nbrimages;
 i=i+1;
 iplus1=mod(i,nbrimages)+1;
 iplus2=mod(i+1,nbrimages)+1;
 [pi,li,ci]=getcommonfeatures(imseq([i,iplus1,iplus2]));
 nbrpoints=length(pi);
 nbrlines=length(li);
 nbrconics=length(ci);
 nbrpc=nbrpoints+nbrconics;

 %check if same center as last time
 if i>1,
  if size(pi,2)==size(pilast,2) & size(ci,2)==size(cilast,2),
   if sum(abs(pi-pilast))+sum(abs(ci-cilast))==0,
    samecenter(i)=1;
   end
  end
 else
  pifirst=pi;
  cifirst=ci;
 end

  % Check for minimality conditions (NB: This is not done strictly)
 if nbrpc<5,
  if nbrpc<3 | nbrlines < 3;
    if i<(nbrimages-1) %closed sequence????
      error(['Too few visible features: ', ...
             num2str(nbrpoints),' points ', ...
             num2str(nbrlines),' lines ', ...
             num2str(nbrconics),' conics']);
    else
      i=nbrimages+1; % no closed sequence
    end
  end
 end
 if i<=nbrimages;
  % get points and centers of conics and lines
  if nbrpoints==0,
    tmp=getconics(imseq{i},ci);
    pts1=pflat(trans{i}*tmp(4:6,:));
    tmp=getconics(imseq{iplus1},ci);
    pts2=pflat(trans{iplus1}*tmp(4:6,:));
    tmp=getconics(imseq{iplus2},ci);
    pts3=pflat(trans{iplus2}*tmp(4:6,:));
  else
    pts1=pflat(trans{i}*getpoints(imseq{i},pi));
    pts2=pflat(trans{iplus1}*getpoints(imseq{iplus1},pi));
    pts3=pflat(trans{iplus2}*getpoints(imseq{iplus2},pi));
    if nbrconics>0,
      tmp=getconics(imseq{i},ci);
      pts1=[pts1,pflat(trans{i}*tmp(4:6,:))];
      tmp=getconics(imseq{iplus1},ci);
      pts2=[pts2,pflat(trans{iplus1}*tmp(4:6,:))];
      tmp=getconics(imseq{iplus2},ci);
      pts3=[pts3,pflat(trans{iplus2}*tmp(4:6,:))];
    end
  end
  lines1=psphere(inv(trans{i})'*gethomogeneouslines(imseq{i},li));
  lines2=psphere(inv(trans{iplus1})'*gethomogeneouslines(imseq{iplus1},li));
  lines3=psphere(inv(trans{iplus2})'*gethomogeneouslines(imseq{iplus2},li));
  [t,c1,c2,c3]=pointslines2tri(pts1,pts2,pts3,lines1,lines2,lines3);
  centers=[centers,[c1(1:2);c2(1:2);c3(1:2)]];
  T=[T;t];
 end %END of if i<=nbrimages
 pilast=pi;
 cilast=ci;
end %END of while

% compare first and last centers
if size(pi,2)==size(pifirst,2) & size(ci,2)==size(cifirst,2),
 if sum(abs(pi-pifirst))+sum(abs(ci-cifirst))==0,
  samecenter(1)=1;
 end
end


P=tri2p(T,nbrimages);

A=cell(1,nbrimages);
for i=1:nbrimages;
  index=2*i-1;
  A{i}=P(index:index+1,:);
end

%calculate last column in camera matrices
M=[];
cindex=2+2*size(centers,2);
sq=1;
changed=0;
for i=[2:size(centers,2),1];
  if samecenter(i)==0 & i~=2,
    cindex=cindex+3;
    changed=1;
  end
  iplus1=mod(i,nbrimages)+1;
  iplus2=mod(i+1,nbrimages)+1;
  M(sq:sq+5,1)=-centers(:,i);
  M(sq:sq+1,2*i:2*i+1)=eye(2);
  M(sq+2:sq+3,2*iplus1:2*iplus1+1)=eye(2);
  M(sq+4:sq+5,2*iplus2:2*iplus2+1)=eye(2);
  if changed==1,
    M(sq:sq+1,cindex:cindex+2)=A{i};
    M(sq+2:sq+3,cindex:cindex+2)=A{iplus1};
    M(sq+4:sq+5,cindex:cindex+2)=A{iplus2};
  end

  sq=sq+6;
end
[u,ss,v]=svd(M);
%tmp=abs(v(1,cindex-1:cindex+2));
%bindex=cindex-2+find(max(tmp)==tmp);
%b=v(:,bindex)/v(1,bindex);

b=v(:,size(v,2));b=b/b(1);

m=motion;
for i=1:nbrimages;
  P=[A{i},b(2*i:2*i+1);0 0 0 1];
  m=addcameras(m,inv(trans{i})*P);
end

s=intseclinear(imseq,m,option);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of main-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,c1,c2,c3]=pointslines2tri(I1,I2,I3,L1,L2,L3);

% [T,c1,c2,c3]=pointslines2tri(I1,I2,I3,L1,L2,L3);
% Calculates the third order affine quasi-tensor
% from point and line correspondences
c1=mean(I1')';
c2=mean(I2')';
c3=mean(I3')';

nbrpoints=size(I1,2);
I1=I1-c1*ones(1,nbrpoints);
I2=I2-c2*ones(1,nbrpoints);
I3=I3-c3*ones(1,nbrpoints);
nbrlines=size(L1,2);

nofconstr=15*nbrpoints+nbrlines;
M=zeros(nofconstr,20);
pointnr=1;
constrnr=1;
while pointnr<=nbrpoints
  M(constrnr,[11,5,2,1])=...
  [I1(1,pointnr) -I1(2,pointnr) I2(1,pointnr) -I2(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[16,10,4,3])=...
  [I1(1,pointnr) -I1(2,pointnr) I3(1,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[20,19,18,17])=...
  [I2(1,pointnr) -I2(2,pointnr) I3(1,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[12,6,3,1])=...
  [I1(1,pointnr) -I1(2,pointnr) I2(1,pointnr) -I3(1,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[13,7,4,1])=...
  [I1(1,pointnr) -I1(2,pointnr) I2(1,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[14,8,3,2])=...
  [I1(1,pointnr) -I1(2,pointnr) I2(2,pointnr) -I3(1,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[15,9,4,2])=...
  [I1(1,pointnr) -I1(2,pointnr) I2(2,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[17,8,6,5])=...
  [I1(1,pointnr) -I2(1,pointnr) I2(2,pointnr) -I3(1,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[18,9,7,5])=...
  [I1(1,pointnr) -I2(1,pointnr) I2(2,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[19,10,7,6])=...
  [I1(1,pointnr) -I2(1,pointnr) I3(1,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[20,10,9,8])=...
  [I1(1,pointnr) -I2(2,pointnr) I3(1,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[17,14,12,11])=...
  [I1(2,pointnr) -I2(1,pointnr) I2(2,pointnr) -I3(1,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[18,15,13,11])=...
  [I1(2,pointnr) -I2(1,pointnr) I2(2,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[19,16,13,12])=...
  [I1(2,pointnr) -I2(1,pointnr) I3(1,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  M(constrnr,[20,16,15,14])=...
  [I1(2,pointnr) -I2(2,pointnr) I3(1,pointnr) -I3(2,pointnr)];
  constrnr=constrnr+1;
  pointnr=pointnr+1;
end
linenr=1;
while linenr<=nbrlines
%  M(constrnr,[6,7,8,9,12,13,14,15])=...
%  [L1(2,linenr)*L2(2,linenr)*L3(2,linenr)...
%   -L1(2,linenr)*L2(2,linenr)*L3(1,linenr)...
%   -L1(2,linenr)*L2(1,linenr)*L3(2,linenr)...
%   L1(2,linenr)*L2(1,linenr)*L3(1,linenr)...
%   -L1(1,linenr)*L2(2,linenr)*L3(2,linenr)...
%   L1(1,linenr)*L2(2,linenr)*L3(1,linenr)...
%   L1(1,linenr)*L2(1,linenr)*L3(2,linenr)...
%   -L1(1,linenr)*L2(1,linenr)*L3(1,linenr)];
  M(constrnr,[6,7,8,9,12,13,14,15])=...
  [L1(1,linenr)*L2(1,linenr)*L3(1,linenr)...
   L1(1,linenr)*L2(1,linenr)*L3(2,linenr)...
   L1(1,linenr)*L2(2,linenr)*L3(1,linenr)...
   L1(1,linenr)*L2(2,linenr)*L3(2,linenr)...
   L1(2,linenr)*L2(1,linenr)*L3(1,linenr)...
   L1(2,linenr)*L2(1,linenr)*L3(2,linenr)...
   L1(2,linenr)*L2(2,linenr)*L3(1,linenr)...
   L1(2,linenr)*L2(2,linenr)*L3(2,linenr)];
  constrnr=constrnr+1;
  linenr=linenr+1;
end
[U,S,V]=svd(M);
[m,n]=size(V);
T=V(:,n)';

function P=tri2p(T,nbrimages);

% P=tri2p(T,nbrimages);
% calculate camera matrices from the third order affine quasi-tensors
% between views i,i+1 and i+2 using the affine closure constraints

noftri=12*size(T,1);
M=zeros(noftri,2*nbrimages);
imnr=1;
trinr=1;
while trinr<noftri
  s=2*imnr-2;

  s3=mod(s+2,2*nbrimages)+1;
  s4=mod(s+3,2*nbrimages)+1;
  s5=mod(s+4,2*nbrimages)+1;
  s6=mod(s+5,2*nbrimages)+1;

  M(trinr,[s+1,s+2,s3,s5])=[T(imnr,12) -T(imnr,6) T(imnr,3) -T(imnr,1)];
  trinr=trinr+1;
  M(trinr,[s+1,s+2,s3,s6])=[T(imnr,13) -T(imnr,7) T(imnr,4) -T(imnr,1)];
  trinr=trinr+1;
  M(trinr,[s+1,s+2,s4,s5])=[T(imnr,14) -T(imnr,8) T(imnr,3) -T(imnr,2)];
  trinr=trinr+1;
  M(trinr,[s+1,s+2,s4,s6])=[T(imnr,15) -T(imnr,9) T(imnr,4) -T(imnr,2)];
  trinr=trinr+1;
  M(trinr,[s+1,s3,s4,s5])=[T(imnr,17) -T(imnr,8) T(imnr,6) -T(imnr,5)];
  trinr=trinr+1;
  M(trinr,[s+1,s3,s4,s6])=[T(imnr,18) -T(imnr,9) T(imnr,7) -T(imnr,5)];
  trinr=trinr+1;
  M(trinr,[s+1,s3,s5,s6])=[T(imnr,19) -T(imnr,10) T(imnr,7) -T(imnr,6)];
  trinr=trinr+1;
  M(trinr,[s+1,s4,s5,s6])=[T(imnr,20) -T(imnr,10) T(imnr,9) -T(imnr,8)];
  trinr=trinr+1;
M(trinr,[s+2,s3,s4,s5])=[T(imnr,17) -T(imnr,14) T(imnr,12) -T(imnr,11)];
  trinr=trinr+1;
M(trinr,[s+2,s3,s4,s6])=[T(imnr,18) -T(imnr,15) T(imnr,13) -T(imnr,11)];
  trinr=trinr+1;
M(trinr,[s+2,s3,s5,s6])=[T(imnr,19) -T(imnr,16) T(imnr,13) -T(imnr,12)];
  trinr=trinr+1;
M(trinr,[s+2,s4,s5,s6])=[T(imnr,20) -T(imnr,16) T(imnr,15) -T(imnr,14)];
  trinr=trinr+1;
  imnr=imnr+1;
end;
[U,S,V]=svd(M);
[m,n]=size(V);
P=V(:,[n-2:n]);

