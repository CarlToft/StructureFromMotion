function plotquadric(quadrics)
%PLOTQUADRIC plotquadric(quadrics) plot quadrics
% INPUT:
%   quadrics - 10xn matrix where columns are dual quadrics

n=20;
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
x = cosphi*cos(theta);
[mx,nx]=size(x); x=x(:)';
y = cosphi*sintheta;
[my,ny]=size(y); y=y(:)';
z = sin(phi)*ones(1,n+1);
[mz,nz]=size(z); z=z(:)';
  X=[x;y;z;ones(size(z))];

for t = 1:size(quadrics,2)
  v = quadrics(:,t);

  l = v2m(v);
  [u,d]=schur(l);
  d=diag(d);
  neg=find(d<-eps);
  if length(neg)>1,
    neg=find(d>eps);
  end

  if length(neg)==1,
    dnew = d;
    dnew(neg)=d(4);
    dnew(4)=d(neg);

    e=eye(4);
    p=e;
    p(:,neg)=e(:,4);
    p(:,4)=e(:,neg);

    d0=sqrt(abs(dnew(1:3)/dnew(4)));
    q=diag([d0;1]);
    T=u*inv(p')*q;
    Xnew=T*X;
    xnew = reshape(Xnew(1,:)./Xnew(4,:),mx,nx);
    ynew = reshape(Xnew(2,:)./Xnew(4,:),my,ny);
    znew = reshape(Xnew(3,:)./Xnew(4,:),mz,nz);
    mesh(xnew,ynew,znew,ones(size(znew)));
  end

end









