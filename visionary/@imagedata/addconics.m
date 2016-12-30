function i=addconics(i,conics,covs)
%IMAGEDATA/ADDCONICS i=addconics(i,conics,covs) adds dual conics
%  and covariances covs

i.conics = [i.conics,conics];
if nargin>=3,
  if iscell(covs)
    i.coniccov = [i.coniccov,covs];
  else
    i.coniccov = [i.coniccov,{covs}];
  end
end




