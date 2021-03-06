function [J,res] = setupLinSystemPosPrior(J,res,P,C_p,numpts,w)
%Bas f�r tangentplanet till rotationsm�ngfalden.
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];


daC = cell(size(P));
dbC = cell(size(P));
dcC = cell(size(P));
dt1C = cell(size(P));
dt2C = cell(size(P));
dt3C = cell(size(P));


nPosConditions = 0;

for i=1:length(P)
    %a,b,c - rotations parametrar f�r kamera i
    %t1,t2,t3 - translations parametrar kameran.
    R0 = P{i}(:,1:3);
    
    
    if all(isfinite(C_p(:,i)))
        daC{i} = Ba*R0*C_p(:,i);
        dbC{i} = Bb*R0*C_p(:,i);
        dcC{i} = Bc*R0*C_p(:,i);
        dtC = eye(3);
        dt1C{i} = dtC(:,1);
        dt2C{i} = dtC(:,2);
        dt3C{i} = dtC(:,3);
        nPosConditions = nPosConditions + 1;
    end
end

if nPosConditions < 2
    error('BA:constraint',  'Not enough camera position priors');
end

if(size(res,1) ~= size(J,1))
    error('BA:constraint', ...
    'Jacobian and residual must have same number of rows');
end

nVars       = 3*numpts + 6*length(P);
iConstraint = 0;
nnzElems    = 3*6*nPosConditions;
nConditions = 3*nPosConditions;
row         = zeros(nnzElems,1);
col         = zeros(nnzElems,1);
data        = zeros(nnzElems,1);
lastentry   = 0;
resAdd      = zeros(nConditions,1);

for i = 1:length(P)

    if all(isfinite(C_p(:,i)))
        row(lastentry+(1:3*6)) = iConstraint+repmat((1:3)',6,1);
        col(lastentry+(1:3*6)) = 3*numpts+(i-1)*6+reshape(repmat(1:6,3,1),18,1);
        data(lastentry+(1:3*6)) = [daC{i}; dbC{i}; dcC{i}; dt1C{i}; dt2C{i}; dt3C{i}];
        lastentry = lastentry+3*6;

        e = P{i}*[C_p(:,i);1];
        resAdd(iConstraint+(1:3)) = e;
        iConstraint = iConstraint+3;
    end

end

J = [J; w*sparse(row,col,data,nConditions, nVars)];
res = [res; w*resAdd];

end