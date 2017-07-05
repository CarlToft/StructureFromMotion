function [J,res] = setupLinSystemRotPrior(J,res,P,R_p,numpts,w)
%Bas f�r tangentplanet till rotationsm�ngfalden.
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];


daR = cell(size(P));
dbR = cell(size(P));
dcR = cell(size(P));

nRotConditions = 0;

for i=1:length(P)
    %ber�kna derivator f�r b�da residualeran i alla bilder
    %a,b,c - rotations parametrar f�r kamera i
    %U1,U2,U3 - 3d punkt parametrar
    %t1,t2,t3 - translations parametrar kameran.
    R0 = P{i}(:,1:3);
    
    
    if ~isempty(R_p{i})
        daR{i} = Ba*R0;
        dbR{i} = Bb*R0;
        dcR{i} = Bc*R0;
        nRotConditions = nRotConditions + 1;
    end
end

if nRotConditions < 2
    error('BA:constraint',  'Not enough camera rotation priors');
end

if(size(res,1)~=size(J,1))
    error('BA:constraint', ...
    'Jacobian and residual must have same number of rows');
end

nVars       = 3*numpts + 6*length(P);
iConstraint = 0;
nnzElems    = 3*9*nRotConditions;
nConditions = 9*nRotConditions;
row         = zeros(nnzElems,1);
col         = zeros(nnzElems,1);
data        = zeros(nnzElems,1);
lastentry   = 0;
resAdd      = zeros(nConditions,1);

for i = 1:length(P)
    if ~isempty(R_p{i})
        row(lastentry+(1:9*3)) = iConstraint+repmat((1:9)',3,1);
        col(lastentry+(1:9*3)) = 3*numpts+(i-1)*6+reshape(repmat(1:3,9,1),27,1);
        data(lastentry+(1:9*3)) = [daR{i}(:); dbR{i}(:); dcR{i}(:)];
        lastentry = lastentry+9*3;
        
        R0 = P{i}(1:3,1:3);
        e = R0-R_p{i};
        resAdd(iConstraint+(1:9)) = e(:);
        iConstraint = iConstraint+9;
    end

end

J = [J; w*sparse(row,col,data,nConditions, nVars)];
res = [res; w*resAdd];

end