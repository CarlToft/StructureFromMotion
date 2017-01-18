function [pos2,residual]=findpatch2(im1,im2,pos1,predpos,precalc)
% [pos2,residual]=findpatch2(im1,im2,pos1,predpos,precalc)
% Finds a patch in image im2,
% using a patch around pos1 in image im1
% by means of affine subpixel normalized correlation
% Input:
%  im1,im2 - image matrices
%  pos1 - position in image im1
%  predpos - predicted position in im2 (optional)
%  precalc - Precalculated data (optional)
%	     This used if findpatch2 is called several times.
% Output:
%  pos2 - position in image im2
%  residual - intensity residual between patches
% See also: findpatch

patchsize=5;

searchsize=120;
stdev=1.5;


if nargin<4
   predpos=pos1;
end;
if sum(pos1(1:2)<patchsize+2)>0 | ...
	pos1(1)>size(im1,2)-patchsize-1 | ...
	pos1(2)>size(im1,1)-patchsize-1,
   residual=Inf;
   pos2=predpos;
else
%%%%%%%%%%%%%%%%%%%
n=2*patchsize+1;
index1=max(1,round(predpos(2))-searchsize):min(size(im2,1),round(predpos(2))+searchsize);
index2=max(1,round(predpos(1))-searchsize):min(size(im2,2),round(predpos(1))+searchsize);
searchpatch=im2(index1,index2);

if nargin<5
  searchones=conv2(searchpatch,ones(n),'valid');
  search2ones=conv2(searchpatch.^2,ones(n),'valid');
  searchsquared=search2ones-searchones.^2/n^2;
  searchsquared(find(searchsquared==0))=1;
else
  searchsquared=precalc(index1(1+patchsize:end-patchsize)-patchsize,...
		 index2(1+patchsize:end-patchsize)-patchsize);

%im2conv=conv2(images{2}.^2,ones(n),'valid')-conv2(images{2},ones(n),'valid').^2/n^2;

end


patch=im1(round(pos1(2))+(-patchsize:patchsize),round(pos1(1))+(-patchsize:patchsize));
patch = flipud(fliplr(patch));

patch=patch-mean(patch(:));
patch2=sum(patch(:).^2);if patch2==0,patch2=1;end;

s_patch=conv2(searchpatch,patch,'valid');

distmap=(searchsquared+patch2-2*s_patch)./sqrt(searchsquared)/sqrt(patch2);

[posxx,posyy]=find(distmap==min(distmap(:)));
posxx=posxx(1)-1+index1(1)+patchsize;
posyy=posyy(1)-1+index2(1)+patchsize;


%pos2=[posyy,posxx]';
%residual=min(distmap(:));

startpos=pos1(1:2)';
destpos=[posyy,posxx];

Astart = calcAstartC(startpos,destpos);

[u, Afirst1,residual] = affineoptC(im1, im2, startpos, Astart, patchsize, stdev);
%[u, Afirst1,residual] = affineoptC(im1, im2, startpos, Astart, 4,2.5);

%figure(1);clf;hold on;plot(imagedata(im1,startpos'),'ro');zoom on;
%figure(2);clf;hold on;plot(imagedata(im2,destpos'),'ro');zoom on;
%residual,
%min(distmap(:))
%keyboard;


um = mean(u');
A00 = [1  0 -um(1) ; ...
       0  1 -um(2) ; ...
       0  0   1    ];
u1 = A00*u;
u1m = std(u1');
A01 = [ 1/u1m(1) 0 0 ; 0 1/u1m(2) 0 ; 0 0 1];
A0 = A01*A00;
Aopt= inv(A0)*Afirst1*A0;

pos2=Aopt*[startpos,1]';





%%%%%%%%%%%%%%%%%%%
end % if pos outside image

%%%%%%%%%%%%%%%%%%%
function [patch, A,f0] = affineoptC(bild1, bild2, pos1, A, patchsize, a);
patch = createpatch(pos1,patchsize);

[f0, Idiff, Idiffgrad] = affineresC(bild1, bild2,patch,A,a);
dA= inf;
counter=0;
while log(norm(dA)+eps)/log(10)>-6 & counter<30;
%  log(norm(dA)+eps)/log(10)
  S = Idiffgrad'*Idiffgrad;
  dA=-inv(S+0.001*eye(size(S,1)))*Idiffgrad'*Idiff;
  Anew=A+reshape([dA ; 0 ; 0 ; 0],size(A))';
  [f1, Idiff, Idiffgrad] = affineresC(bild1, bild2,patch,Anew,a);
  while f1>f0 & counter<30,
%    f0
%    f1
%    disp('Det har tar vi en gang till');
    dA=dA/2;
    Anew=A+reshape([dA ; 0 ; 0 ; 0],size(A))';
    [f1, Idiff, Idiffgrad] = affineresC(bild1, bild2,patch,Anew,a);
    counter=counter+1;
  end;

  if f1<=f0,
    f0=f1;
    A = Anew;
  end;
  counter=counter+1;
end;


function [f,intendiff_normalized, intendiffgrad, i1, i2]=...
		affineresC(image1,image2,points1,A,a);

% First calculate A0

um = mean(points1');

A00 = [1  0 -um(1) ; ...
       0  1 -um(2) ; ...
       0  0   1    ];

u1 = A00*points1;
u1m = std(u1');
A01 = [ 1/u1m(1) 0 0 ; 0 1/u1m(2) 0 ; 0 0 1];
A0 = A01*A00;

points2 = inv(A0)*A*A0*points1;

if max(points2(1,:))>(size(image2,2)-5) | ...
   max(points2(2,:))>(size(image2,1)-5) | ...
   min(points2(1,:))<6 | ...
   min(points2(2,:))<6,

   f=Inf;
   intendiff_normalized=Inf;
   intendiffgrad=Inf*ones(1,6);
   i1=Inf;
   i2=Inf;
else

%%%%%%%



[i1, dummy] = measure(image1,points1,a);
[i2, i2grad] = measure(image2,points2,a);
i1mean=i1-mean(i1);
i2mean=i2-mean(i2);

intendiff = (i2mean-i1mean);
i1square=sum(i1mean.^2);
N=i1square*sum(i2mean.^2);

M = eye([9 9]);
d1 = inv(A0)*reshape(M(1,1:9),3,3)'*A0;
d2 = inv(A0)*reshape(M(2,1:9),3,3)'*A0;
d3 = inv(A0)*reshape(M(3,1:9),3,3)'*A0;
d4 = inv(A0)*reshape(M(4,1:9),3,3)'*A0;
d5 = inv(A0)*reshape(M(5,1:9),3,3)'*A0;
d6 = inv(A0)*reshape(M(6,1:9),3,3)'*A0;

dAx1 = d1 * points1;
dAx2 = d2 * points1;
dAx3 = d3 * points1;
dAx4 = d4 * points1;
dAx5 = d5 * points1;
dAx6 = d6 * points1;

tmp = sum(i2grad .* dAx1(1:2,:));
dfel1 = tmp-mean(tmp);
tmp = sum(i2grad .* dAx2(1:2,:));
dfel2 = tmp-mean(tmp);
tmp = sum(i2grad .* dAx3(1:2,:));
dfel3 = tmp-mean(tmp);
tmp = sum(i2grad .* dAx4(1:2,:)); 
dfel4 = tmp-mean(tmp);
tmp = sum(i2grad .* dAx5(1:2,:));
dfel5 = tmp-mean(tmp);
tmp= sum(i2grad .* dAx6(1:2,:));
dfel6 = tmp-mean(tmp);


dN1=2*sum(i2mean.*dfel1)*i1square;
ddfel1=dfel1/N^(1/4)-intendiff/4/N^(5/4)*dN1;
dN2=2*sum(i2mean.*dfel2)*i1square;
ddfel2=dfel2/N^(1/4)-intendiff/4/N^(5/4)*dN2;
dN3=2*sum(i2mean.*dfel3)*i1square;
ddfel3=dfel3/N^(1/4)-intendiff/4/N^(5/4)*dN3;
dN4=2*sum(i2mean.*dfel4)*i1square;
ddfel4=dfel4/N^(1/4)-intendiff/4/N^(5/4)*dN4;
dN5=2*sum(i2mean.*dfel5)*i1square;
ddfel5=dfel5/N^(1/4)-intendiff/4/N^(5/4)*dN5;
dN6=2*sum(i2mean.*dfel6)*i1square;
ddfel6=dfel6/N^(1/4)-intendiff/4/N^(5/4)*dN6;


intendiffgrad = [ ddfel1 ; ddfel2 ; ddfel3 ; ddfel4 ; ddfel5 ; ddfel6]';

intendiff_normalized=intendiff'/N^(1/4);


f = intendiff_normalized'*intendiff_normalized;

end;

function u = createpatch(pos,s);
% Creates a list of points that make up a patch.
% (Rectangular shape)

i=1;
for x=-s:s
  for y=-s:s
    u(1:3,i)=[pos(1)+x pos(2)+y 1]';
    i=i+1;
  end;
end;

function [inten, intengrad]=measure(image,points,a);

nbr=size(points,2);

inten=zeros(1,nbr);
intengrad=zeros(2,nbr);

[m,n]=size(image);
NN=round(a*3);


for i=1:nbr,


%%%%%%%%  [In, Igrad] = measurepoint(image,points(1:2,i),a);
%measurepoint
%function [inten, intengrad]=measurepoint(bild,punkt,a);


x0=points(1,i);
y0=points(2,i);


cutx=max(round(x0)-NN,1):min(round(x0)+NN,n);
cuty=max(round(y0)-NN,1):min(round(y0)+NN,m);
cutbild=image(cuty,cutx);

[x,y]=meshgrid(cutx,cuty);

filter=exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^2*pi);
filtx=2*(x-x0).*exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^4*pi);
filty=2*(y-y0).*exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^4*pi);
%if size(filtx,1)~=size(cutbild,1) | size(filtx,2)~=size(cutbild,2),
%  keyboard;
%end
Igrad(1,1)=sum(sum(filtx.*cutbild));
Igrad(2,1)=sum(sum(filty.*cutbild));
In=sum(sum(filter.*cutbild));

%measurepoint
%%%%%%%%%%%%%%%%%%%%


  inten(i) = In;
  intengrad(1:2,i) = Igrad;
end;

function Astart = calcAstartC(pos1, pos2);

delta = pos2-pos1;
A = eye(3);
A(1,3)=delta(1);
A(2,3)=delta(2);

um=pos1;

A00 = [1  0 -um(1) ; 0  1 -um(2) ;  0  0   1    ];
A01 = [ 1/3.1754 0 0 ; 0 1/3.1754 0 ; 0 0 1];
A0 = A01*A00;

%Astart = inv(A0)*A*A0;
Astart = A0*A*inv(A0);


