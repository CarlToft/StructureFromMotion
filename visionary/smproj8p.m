function [s,m]=smproj8p(imseq,index);
% [s,m]=smproj8p(imseq,index) solves the projective structure and motion
% problem for image-point correspondences.
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

% normalize coordinates for numerical conditioning
H1 = getnormtrans(imseq{1});
H2 = getnormtrans(imseq{2});

p0=getcommonfeatures(imseq(index));
nbrpoints=length(p0);

if nbrpoints<8,
  error('Too few points visible in all images');
end

x1 = getpoints(imseq{index(1)},p0);
x2 = getpoints(imseq{index(2)},p0);
u1 = H1*x1;
u2 = H2*x2;
[m,s]=sm3(u1,u2);
F = ptobi(m);
F = H1 * F * H2; % undo coordinate normalization
m = bitop(F);

function [m,s]=sm3(u1,u2);
% [P1,P2,U,l1,l2,u1t,u2t]=sm3(u1,u2);
% Calculate structure and motion given
% at least 8 corresponding image points in two images
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
% 3x3 fundamental Matrix F

% We are now going to construct a matrix ev
% containing these linear constraints on the
% parameters of the fundamental matrix F
% such that ev*F(:)=0.
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
F=reshape(u(:,9),3,3);

% Make sure that F has rank 2
[u s v] = svd(F);
s(end) = 0;
F = u*s*v';

% Also renormalize so that F has unit norm in the frobenius.
F=F/norm(F,'fro');

% compute motion and structure objects
m = bitop(F);
s = intsec2views(m, {imagedata([],u1), imagedata([],u2)});