function i=addpoints(i,points,covs)
%IMAGEDATA/ADDPOINTS i=addpoints(i,points,cov) adds homog. or cart. points
%  and covariances covs

if size(points,1) == 2,
  i.points = [i.points,[points;ones(1,size(points,2))]];
else
  i.points = [i.points,points];
end

if nargin>=3;
  if iscell(covs),
    i.pointcov=[i.pointcov,covs];
  else
    i.pointcov=[i.pointcov,{covs}];
  end
end

