function [dpointvar, dcamvar] = restoreFixedVars(...
    d,nPoints,nCams, fixedCams, fixedPoint, fixAllPoints)

    firstCamI = 3*nPoints+1;
    if fixAllPoints
        dpointvar = zeros(3*nPoints,1);
        firstCamI = 1;
    elseif fixedPoint
        dpointvar = [0; d(1:(3*nPoints-1))];
        firstCamI = firstCamI-1;
    else
        dpointvar = d(1:(3*nPoints));
    end

    if isempty(fixedCams)
        dcamvar = d(firstCamI:end);
    else
        dcamvar = zeros(nCams*6,1);
        index = setdiff(1:nCams,fixedCams);
        keepcols = [(index-1)*6+1;(index-1)*6+2;(index-1)*6+3;...
                    (index-1)*6+4;(index-1)*6+5;(index-1)*6+6];
        keepcols = keepcols(:)';
        dcamvar(keepcols) = d(firstCamI:end);
    end

    dpointvar = reshape(dpointvar, [3 nPoints]);
    dcamvar = reshape(dcamvar,[6 nCams]);
end