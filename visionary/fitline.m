function i=fitline(points,normals);
% i=fitline(points,normals);
% INPUT:
%   points - image point positions in homogeneous coordinates (3xn matrix).
%   normals - error normals for each point (3xn matrix)
% OUTPUT:
%   i - imagedata object
% If no normals are specified, the error is assumed to be isotropic.


if nargin<=1,
  [l,L]=fitline1(points);
else
  [l,L]=fitline2(points,normals);
end

i=imagedata;
i=addlines(i,l,L);


function [l,L]=fitline1(points);
% [l,L]=fitline1(points);
% INPUT:
%   points - image point positions in homogeneous coordinates (3xn matrix).
% OUTPUT:
%   l - line parameters
%   L - cholesky factorisation of inverse of covariance of l.
%

% Normalise points
  points=pflat(points);

% Rough estimate
  [u,S,v]=svd(points);
  l=u(:,3);
  l=l/norm(l(1:2));

% Find line such that l'*idealpoints=0, with idealpoints as close
% to points as possible.
  for i=1:5;
   dldx = [0 0 100;-l(2) l(1) 0]';
%   res = (points'*l)./stddevs';
%   dresdx = points'./(stddevs'*ones(1,3))*dldx;
   res = (points'*l);
   dresdx = points'*dldx;
   dl=-dldx*(dresdx\res);
   l=l+dl;
   l=l/norm(l(1:2));
%   norm(dl)
  end;
  estimstd = sqrt((res'*res) / (size(res,1)-2));

  Fx = dresdx'*dresdx/estimstd^2;
  Cx = pinv(Fx);
%  Lx = chol(Fx);
%  Fl = dldx*Fx*dldx'; % FEL
%  L = Lx*dldx';
  Cl = dldx*Cx*dldx';
  %pinv(Fl);
  [U,S,V]=svd(Cl);
  U=U(:,1:2);
  S=S(1:2,1:2);
  L=inv(sqrtm(S))*U';

function [l,L]=fitline2(points,normals);
% [l,L]=fitline2(points);
% INPUT:
%   points - image point positions in homogeneous coordinates (3xn matrix).
%   normals - error normals for each point (3xn matrix)
% OUTPUT:
%   l - line parameters
%   L - cholesky factorisation of inverse of covariance of l.
%

  [u,S,v]=svd(points);
  u=u(:,3);
  for k=1:5;
   u=u/norm(u);
   ss = -(u'*points)./(u'*normals);
   normss=norm(ss);
   idealpoints=points+normals*diag(ss);
   meanpoint = mean(points')';
   dssdu = -(points*diag(ones(size(ss))./(u'*normals)) - ...
   normals*diag(u'*points./((u'*normals).^2)));
   [U,S,V]=svd(dssdu');
   du = V(:,1:2)*inv(S(1:2,1:2))*U(:,1:2)'*(ss');
   u=u-du;
  end;
%  normdu=norm(du)
  stddev=sqrt(ss*ss'/(size(ss,2)-2));
  M= V(:,1:2)*inv(S(1:2,1:2))*U(:,1:2)';
  Cu = M*diag(stddev^2*ones(size(ss)))*M';
%  [U,S,V]=svd(Cu);
%  U=U(:,1:2);
%  S=S(1:2,1:2);
%  L=inv(sqrtm(S))*U';
%  keyboard;

  l=u/norm(u(1:2));
  dldu = eye(3)/norm(u(1:2)) -u*(1/norm(u(1:2))^3)*[u(1:2)' 0];
  Cl = dldu*Cu*dldu';
  [U,S,V]=svd(Cl);
  U=U(:,1:2);
  S=S(1:2,1:2);
  L=inv(sqrtm(S))*U';
