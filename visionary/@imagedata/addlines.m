function i=addlines(i,lines,covs)
%IMAGEDATA/ADDLINES i=addlines(i,lines,covs) adds homog. lines
%  and covariances covs

if size(i.lines,2)==0,
  i.lines = lines;
else
  i.lines = [i.lines,lines];
end

if nargin>=3,
  if iscell(covs)
    i.linecov = [i.linecov,covs];
  else
    i.linecov = [i.linecov,{covs}];
  end
end


