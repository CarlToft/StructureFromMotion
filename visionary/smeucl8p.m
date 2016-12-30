function [s,m]=smeucl8p(imseq,index);
% [s,m]=smeucl8p(imseq,index) solves the structure and motion
% problem for points
% INPUT:
%   imseq - cell list containing imagedata objects
%   index - reference to 2 view (optional)
% OUTPUT:
%   s - structure
%   m - motion
% NB: At least 8 points must be visible in all images

if nargin<2,
    index=[1,2];
end

p0=getcommonfeatures(imseq(index));
nbrpoints=length(p0);

if nbrpoints<8,
  error('Too few points visible in all images');
end

u1 = getpoints(imseq{index(1)},p0);
u2 = getpoints(imseq{index(2)},p0);

[m,s]=sm3(u1,u2);



function [m,s]=sm3(u1,u2);
% [P1,P2,U,l1,l2,u1t,u2t]=sm3(u1,u2);
% Calculate structure and motion given
% at least 8 corresponding image points in two images
% assuming that the point coordinates have been
% corrected for internal calibration.
% Input: u1, u2 matrices of type 3x6. Each column represents
%        a point. Corresponding columns in u1 and u2 are
%        the image points of the same point in image 1 and 2.
% Output: P1, P2 - 3x4 camera matrices (camera position and orientation)
%         U - 4x6 matrix representing the position of the six points.
%         l1, l2 - depth vectors.
%         u1t, u2t - reprojected points.
% We have u1t*diag(l1) = P1*U and 
% and we hope to have u1t approximately equal to u1,
% and similarily for image 2.
%
% Copyright Fredrik Kahl, Kalle Åström 2006.


% The corresponding points
% give linear constraints on the
% 3x3 essential Matrix E

% We are now going to construct a matrix ev
% containing these linear constraints on the
% parameters of the essential matrix E
% such that ev*E(:)=0.
ev=[];
nr_of_points=size(u1,2);
if nr_of_points<8,
  error('The routine sm3 requires at least 8 points');
end;
for i=1:nr_of_points
 v=u1(:,i)*u2(:,i)';
 ev=[ev v(:)];
end;
[u,s,v]=svd(ev);
% Here one ought to test if the matrix ev has 8 sufficiently
% large singular values and that the last singular is sufficiently
% small
% If so then the last column of u in the singular value 
% decomposition gives the vector uu that minimized norm(ev*u,'fro');
E=reshape(u(:,9),3,3);

% This does not guarantee that E is an essential matrix.
% Therefore we first find a matrix that is close to E but
% fullfills the constraints on being an essential matrix.
% i.e. renormalise so that E has singular values 1, 1 and 0.
[u,s,v]=svd(E);
smean=(s(1,1)+s(2,2))/2;
s=diag([smean smean 0]);
E=u*s*v';
% Also renormalize so that E has unit norm in the frobenius.
E=E/norm(E,'fro');

% Then we optimize so that E is essential while minimizing
% norm(ev*u,'fro');
% This is done by iterative techniques?

% Calculate camera matrices from E.
[t,R]=invessentialm(E,u1,u2);
% There are several possibilities
% P1=[I 0], P2=[R1' (-R1'*t)];
% P1=[I 0], P2=[R1' (+R1'*t)];
% P1=[I 0], P2=[R2' (-R2'*t)];
% P1=[I 0], P2=[R2' (+R2'*t)];

P1=[eye(3) zeros(3,1)];
P2=[R' -R'*t];
m=motion({P1,P2});
s=intsec2views(m,{imagedata([],u1),imagedata([],u2)});

% Calculate the structure U.
%for i=1:nr_of_points,
% M=[P1 u1(:,i) zeros(3,1); P2 zeros(3,1) u2(:,i)];
 % One ought to check that the intersection is well defined.
% [u,s,v]=svd(M);
%diag(s)'
% U(1:4,i)=v(1:4,6);
% l1(1,i)=-v(5,6);
% l2(1,i)=-v(6,6);
%end;
%[U,aa]=pflat(U);
%[u1t,l1]=psphere(P1*U);
%[u2t,l2]=psphere(P2*U);
%o1=sign(sum(u1t(:,1:nr_of_points).*u1(:,1:nr_of_points)));
%o2=sign(sum(u2t(:,1:nr_of_points).*u2(:,1:nr_of_points)));
% Perhaps one ought to check that u1t = u1 and that u2t=u2.
