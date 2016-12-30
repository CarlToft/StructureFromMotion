function [pos2,residual]=affineoptimize(im1, im2, startpos, destpos, patchsize, stdev, Astart,maxsteps);
%local patch optimization over affine transformations
% with some regularity (first 2x2 submatrix should approximately be eye(2))
if nargin<7,
  Astart=eye(2);
end
if nargin<8,
  maxsteps=30;
end
startpos=[startpos(1);startpos(2);1];
destpos=[destpos(1);destpos(2);1];

borderignore=round(3*stdev+2*patchsize+6);
[imsz1,imsz2]=size(im2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate patchpositions
%%%%%%%%%%%%%%%%%%%%%%%%%%%


cnt=1;
for x=-patchsize:patchsize,
  for y=-patchsize:patchsize,
    patchpos(1:3,cnt)=[startpos(1)+x,startpos(2)+y 1]';
    cnt=cnt+1;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate start transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta = destpos-startpos;
A = eye(3);
A(1,3)=delta(1);
A(2,3)=delta(2);
A00 = [1,0,-startpos(1);0,1,-startpos(2);0,0,1];
u1m = std((A00*patchpos)');
A0 = [ 1/u1m(1) 0 0 ; 0 1/u1m(2) 0 ; 0 0 1]*A00;
A=A0*A*inv(A0);
A(1:2,1:2)=Astart;

gscale=1;
detA=det(A(1:2,1:2));
i1 = measure(im1,patchpos,stdev);


[i2, i2grad] = measure(im2,inv(A0)*A*A0*patchpos,stdev);
[f0, Idiff, Idiffgrad] = affineres(i1,i2,i2grad,A,A0,patchpos,detA,gscale);

dA=inf;f1=inf;
counter=0;
while log(norm(dA)+eps)/log(10)>-6 & counter<maxsteps;
%  log(norm(dA)+eps)/log(10)
  S = Idiffgrad'*Idiffgrad;
  dA=-inv(S+0.0001*eye(size(S,1)))*Idiffgrad'*Idiff;
  Anew=A+reshape([dA(1:6) ; 0 ; 0 ; 0],3,3)';gscalenew=gscale+dA(7);

  Aopt=inv(A0)*Anew*A0; pos2=Aopt*startpos;
  if  pos2(1)>=borderignore & pos2(1)<=imsz2-borderignore & ...
      pos2(2)>=borderignore & pos2(2)<=imsz1-borderignore,
    [i2, i2grad] = measure(im2,Aopt*patchpos,stdev);
    [f1, Idiff, Idiffgrad] = affineres(i1,i2,i2grad,Anew,A0,patchpos,detA,gscalenew);
  else
    counter=maxsteps;
  end

  while counter<maxsteps & f1>f0,

    if pos2(1)<borderignore | pos2(1)>imsz2-borderignore | ...
       pos2(2)<borderignore | pos2(2)>imsz1-borderignore,

      counter=maxsteps;
    else
      dA=dA/2;
%    f0,f1,disp('Det har tar vi en gang till');

      Anew=A+reshape([dA(1:6) ; 0 ; 0 ; 0],3,3)';gscalenew=gscale+dA(7);
      Aopt=inv(A0)*Anew*A0; pos2=Aopt*startpos;
      [i2, i2grad] = measure(im2,Aopt*patchpos,stdev);
      [f1, Idiff, Idiffgrad] = affineres(i1,i2,i2grad,Anew,A0,patchpos,detA,gscalenew);
      counter=counter+1;
    end
  end;

  if f1<=f0,
    f0=f1;
    A = Anew;gscale=gscalenew;
  end;
  counter=counter+1;
end;

residual=sqrt(f0)/(2*patchsize+1);


%%%end of affineoptimize


function [f,intendiff, intendiffgrad]=...
		affineres(i1,i2,i2grad,A,A0,patchpos,detA,gscale);

intendiff = (i2*gscale-i1)';

M = eye([9 9]);
d1 = inv(A0)*reshape(M(1,1:9),3,3)'*A0;
d2 = inv(A0)*reshape(M(2,1:9),3,3)'*A0;
d3 = inv(A0)*reshape(M(3,1:9),3,3)'*A0;
d4 = inv(A0)*reshape(M(4,1:9),3,3)'*A0;
d5 = inv(A0)*reshape(M(5,1:9),3,3)'*A0;
d6 = inv(A0)*reshape(M(6,1:9),3,3)'*A0;

dAx1 = d1 * patchpos;
dAx2 = d2 * patchpos;
dAx3 = d3 * patchpos;
dAx4 = d4 * patchpos;
dAx5 = d5 * patchpos;
dAx6 = d6 * patchpos;

dfel1 = sum(i2grad .* dAx1(1:2,:));
dfel2 = sum(i2grad .* dAx2(1:2,:));
dfel3 = sum(i2grad .* dAx3(1:2,:));
dfel4 = sum(i2grad .* dAx4(1:2,:)); 
dfel5 = sum(i2grad .* dAx5(1:2,:));
dfel6 = sum(i2grad .* dAx6(1:2,:));

dfel7 = i2;


intendiffgrad = [ dfel1 ; dfel2 ; dfel3 ; dfel4 ; dfel5 ; dfel6 ; dfel7]';


% add regularization det(A(1:2,1:2))=1

weight=100;
intendiff=[intendiff;weight*[A(1,1)-1;A(1,2);A(2,1);A(2,2)-1;gscale-1]];
intendiffgrad=[intendiffgrad;...
          weight*[1,0,0,0,0,0,0];...
          weight*[0,1,0,0,0,0,0];...
          weight*[0,0,0,1,0,0,0];...
          weight*[0,0,0,0,1,0,0];...
          weight*[0,0,0,0,0,0,1]];

f = intendiff'*intendiff;



function [inten, intengrad]=measure(image,points,a);

nbr=size(points,2);

inten=zeros(1,nbr);
intengrad=zeros(2,nbr);

NN=round(a*3);


for i=1:nbr,


%%%%%%%%  [In, Igrad] = measurepoint(image,points(1:2,i),a);
%measurepoint
%function [inten, intengrad]=measurepoint(bild,punkt,a);


x0=points(1,i);
y0=points(2,i);


cutx=(round(x0)-NN):(round(x0)+NN);
cuty=(round(y0)-NN):(round(y0)+NN);
%cutx=max(round(x0)-NN,1):min(round(x0)+NN,n);
%cuty=max(round(y0)-NN,1):min(round(y0)+NN,m);
cutbild=image(cuty,cutx);

[x,y]=meshgrid(cutx,cuty);

filter=exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^2*pi);
filtx=2*(x-x0).*exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^4*pi);
filty=2*(y-y0).*exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^4*pi);

Igrad(1,1)=sum(sum(filtx.*cutbild));
Igrad(2,1)=sum(sum(filty.*cutbild));
In=sum(sum(filter.*cutbild));

%measurepoint
%%%%%%%%%%%%%%%%%%%%


inten(i) = In;
intengrad(1:2,i) = Igrad;



end






