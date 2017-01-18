function [l,L]=fitline(points,stddevs);
% [l,L]=fitline(points,stddevs);
% INPUT:
%   points - image point positions in homogeneous coordinates (3xn matrix).
%   stddevs - estimated standard deviation for each image point (1xn matrix)
% OUTPUT:
%   l - line parameters
%   L - cholesky factorisation of covariance of l.
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
   res = (points'*l)./stddevs';
   dresdx = points'./(stddevs'*ones(1,3))*dldx;
   dl=-dldx*(dresdx\res);
   l=l+dl;
   l=l/norm(l(1:2));
   norm(dl)
  end;
  Fx = dresdx'*dresdx;
  Lx = chol(Fx)
%  Fl = dldx*Fx*dldx';
  L = Lx*dldx'


