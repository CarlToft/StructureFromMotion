function m=reseclinear(imseq,s,v,option);
% RESECLINEAR m=reseclinear(imseq,s,v,option) resection, linear method using points
% INPUT:
%   imseq - imagedata array
%   s - structure
%   v - views (If not specified all views are resected)
%   option: 'nocoordtrans' - No coordinate transformation
% OUTPUT:
%   m - motion, i.e. resected cameras

if nargin <=2 | isempty(v),
  v=1:length(imseq);
end
if nargin <=3,
  option = '';
elseif ~strcmp('nocoordtrans',option),
  error('Unknown option');
end
m=motion; % create empty motion object
U = psphere(getpoints(s));
good=isfinite(U(1,:));

for i=v;
  if strcmp('nocoordtrans',option),
    T=eye(3);
  else
    T = getnormtrans(imseq{i});
  end
  u = psphere(T*getpoints(imseq{i}));
  good2=isfinite(u(1,:));
  index=find(good & good2);
  p=inv(T)*camerafromuandU(u(:,index),U(:,index));
  p=sign(det(p(:,1:3)))*p/norm(p);
  m=addcameras(m,p);
end;

%end of main function

function p=camerafromuandU(u,U);
% p=camerafromuandU(u,U) - calculates p
% given 3D object points in homogeneous coordinates U
% and 2D image points in homogeneous coordinates u

nbrpoints=size(U,2);

M=[];
j=0;
for i=1:nbrpoints;
  if sum(~isfinite(u(:,i)))==0,
    Up=[U(:,i)',zeros(1,8);zeros(1,4),U(:,i)',zeros(1,4);zeros(1,8),U(:,i)'];
    M=[M,zeros(j*3,1);Up,zeros(3,j),u(:,i)];
    j=j+1;
  end
end

[uq,sq,vq]=svd(M);
w=vq(1:12,12+j);
p = 	[w(1:4)';w(5:8)';w(9:12)'];







