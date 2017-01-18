% Matches two sets of SIFT vectors. Takes as argument two matrices with 128
% rows (SIFT features as columns), and returns two row vectors ind1 and
% ind2, where all pairs (ind1{i}, ind2{i}) correspond to matched vectors. 
function [ind1, ind2] = matchsiftvectors(SIFT1, SIFT2, loweRatio)
    tmp1 = normc(double(SIFT1));
    tmp2 = normc(double(SIFT2));
    
    kdtree = vl_kdtreebuild(tmp2,'NumTrees', 12);
    [index,~] = vl_kdtreequery(kdtree,tmp2,tmp1,'MAXCOMPARISONS',50,'NUMNEIGHBORS',2);
    maxvals = sum(tmp2(:,index(1,:)).*tmp1,1)';
    secondmaxvals = sum(tmp2(:,index(2,:)).*tmp1,1)';
    tmp = acos(maxvals) < loweRatio*acos(secondmaxvals);
    
    matches = zeros(1,size(tmp1,2)); 
    matches(tmp) = index(1,tmp);
    ind1 = find(matches ~= 0);
    ind2 = matches(ind1);
    
    % Remove non unique matches (Caused by NN approximation?)
    [~,I,~]=unique(ind2);
    ind1 = ind1(I);
    ind2 = ind2(I);
end