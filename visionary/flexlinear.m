function T=flexlinear(mot,pp);
% function T=flexlinear(mot,pp);
% Computes projective to Euclidean upgrade linearly
% using all intrinsics known but focallengths
% Input:
%   mot - projective motion
%   pp - approximate principal point (optional)
% Output:
%   T - 4 x 4 coordinate transformation for Euclidean upgrade
%       Use changecsystem(mot,T) for upgrade!
% Default: pp=[0,0]';

if nargin<2,
  pp=[0,0]';
end;

Kpp=eye(3);Kpp(1:2,3)=-pp;

nbrcameras=size(mot);

if nbrcameras<3
  error('Too few cameras');
end;


P=zeros(4,3*nbrcameras);
for i=1:nbrcameras,
  P(:,3*i-2:3*i)=(Kpp*getcameras(mot,i))';
end

[u,s,v]=svd(P);
  
T0=inv(u*s(:,1:4));
P=T0*P;
M=zeros(4*nbrcameras,10);

for i=1:nbrcameras;

  p=P(:,3*i-2:3*i)';
  p=p/norm(p);

  ii=(i-2)*4+1;

  M(4*i-3:4*i,:)=[...
   p(1,1)*p(1,1)-p(2,1)*p(2,1),-2*p(2,1)*p(2,2)+2*p(1,1)*p(1,2),p(1,2)*p(1,2)-p(2,2)*p(2,2),-2*p(2,1)*p(2,3)+2*p(1,1)*p(1,3),2*p(1,2)*p(1,3)-2*p(2,2)*p(2,3),p(1,3)*p(1,3)-p(2,3)*p(2,3),-2*p(2,1)*p(2,4)+2*p(1,1)*p(1,4),2*p(1,2)*p(1,4)-2*p(2,2)*p(2,4),-2*p(2,3)*p(2,4)+2*p(1,3)*p(1,4),p(1,4)*p(1,4)-p(2,4)*p(2,4);
   p(1,1)*p(2,1),p(2,1)*p(1,2)+p(2,2)*p(1,1),p(2,2)*p(1,2),p(2,1)*p(1,3)+p(2,3)*p(1,1),p(2,2)*p(1,3)+p(2,3)*p(1,2),p(2,3)*p(1,3),p(2,1)*p(1,4)+p(2,4)*p(1,1),p(2,2)*p(1,4)+p(2,4)*p(1,2),p(2,3)*p(1,4)+p(2,4)*p(1,3),p(2,4)*p(1,4);
   p(1,1)*p(3,1),p(3,1)*p(1,2)+p(3,2)*p(1,1),p(3,2)*p(1,2),p(3,1)*p(1,3)+p(3,3)*p(1,1),p(3,2)*p(1,3)+p(3,3)*p(1,2),p(3,3)*p(1,3),p(3,1)*p(1,4)+p(3,4)*p(1,1),p(3,2)*p(1,4)+p(3,4)*p(1,2),p(3,3)*p(1,4)+p(3,4)*p(1,3),p(3,4)*p(1,4);
   p(2,1)*p(3,1),p(3,1)*p(2,2)+p(3,2)*p(2,1),p(3,2)*p(2,2),p(3,1)*p(2,3)+p(3,3)*p(2,1),p(3,2)*p(2,3)+p(3,3)*p(2,2),p(3,3)*p(2,3),p(3,1)*p(2,4)+p(3,4)*p(2,1),p(3,2)*p(2,4)+p(3,4)*p(2,2),p(3,3)*p(2,4)+p(3,4)*p(2,3),p(3,4)*p(2,4)];

end





[u,s,v]=svd(M);
Om=v2m(v(:,10));
[u,s,v]=svd(Om);
%dd=diag(s);
%dd(4)=0;
%Om=u*diag(dd)*v';

T=inv(T0'*u*diag(sqrt([s(1,1),s(2,2),s(3,3),1])));
T=T/norm(T);




















