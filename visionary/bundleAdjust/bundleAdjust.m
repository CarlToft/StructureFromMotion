function [P, U, res] = bundleAdjust(P, varargin)

parser = inputParser;
parser.FunctionName = mfilename;
parser.StructExpand = false;

defaultu.pointnr = 0;
defaultu.points = {};
defaultu.index = {};

addRequired(parser, 'P', @isCameraMatrixCellArray);
addOptional(parser, 'U', zeros(4,0), @(x) isreal(x) && 4==size(x,1));
addOptional(parser, 'tracks', defaultu, @(x) isPointTracks(P,x));
addOptional(parser, 'iterations', 20, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addOptional(parser, 'lambda', 0.001, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'relPoses', [], @(x) isRelPoseArray(P,x));
addParameter(parser, 'absPoses', [], @(x) isAbsPoseArray(P,x));
addParameter(parser, 'absPos', [], @(x) isreal(x) && size(x,1) == 3);
addParameter(parser, 'fixedCams', [], @(x) isCamArray(P,x));
addParameter(parser, 'fixedPoint', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'fixAllPoints', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'minLambda', 1e-6, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'maxLambda', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'minResChange', 1e-4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'verbosity', 0, @(x) isaninteger(x,0,2) && isscalar(x));
addParameter(parser, 'schur', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'camWeight', 1/0.1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'posWeight', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);

parse(parser,P,varargin{:});
parm = parser.Results;

%Further parameter checks
assert(size(parm.U,2) == parm.tracks.pointnr, ...
    'Number of 3d points must match between U and tracks');

assert(~parm.fixedPoint || parm.tracks.pointnr > 0, ...
    'There are no points to fix');

assert(~parm.fixAllPoints || parm.tracks.pointnr > 0, ...
    'There are no points to fix');

%Using Schur complement to eliminate points first only makes sense if there
%are any points
parm.schur = parm.schur && parm.tracks.pointnr > 0 && ~parm.fixAllPoints;
useProgressBar = exist('progressBar', 'file') > 0 && parm.verbosity==1;

%There is a scheme for lambda where it decreases from the initial value
%by dividing by 1.5 every iteration
lambda = min(parm.maxLambda,parm.lambda);
nCams = length(P);
nPoints = parm.tracks.pointnr;

for i = 1:nCams
    [KKK,R] = rq(P{i}(:,1:3));
    KK{i} = KKK;
    P{i} = KKK\P{i};
    if nPoints > 0
        parm.tracks.points{i} = pflat(KKK\parm.tracks.points{i});
    end
end
parm.P = P;

Unew = parm.U;
Pnew = parm.P;
if parm.verbosity>=2
    fprintf('#\tResidual\tLambda\n')
end
%Setup and solve linear system and refine the linearization iteratively
for i = 1:parm.iterations
    %Setup a linear system (linearized in the current guess of P and U)
    [B,A] = setupLinSystem(parm);
    
    res = B'*B;
    if i==1 && parm.verbosity>=2
        fprintf('%d\t%f\t%f',0,res,lambda)
    end
    resnew = Inf;
    %If lambda is very small, step size may become too large and go
    %outside the "trust region" so increase lambda until the new residual
    %is smaller than the old.
    while resnew > res
        %If lambda is very large, step size will be very small and almost 
        %no progress towards the minumum will be made.
        if lambda > parm.maxLambda
            if parm.verbosity>=2
                fprintf('\nLambda has grown too large, exiting\n'); 
            end
            %Assign return parameters
            P = parm.P;
            U = parm.U;
            for k = 1:length(P)
                P{k} = KK{k}*P{k};
            end
            return
        end
        if ~parm.schur
            %Marquardt update
            C = A'*A;
            ATA_diag = spdiags(sqrt(abs(diag(C))),0,sparse(size(A,2),size(A,2)));
            C = (C+lambda*ATA_diag);
            c = A'*B;
            if parm.verbosity>=2
                fprintf('\t\tSolving mod-Newton system.\n');
            end
            d = -C\c;
        else
            % Eliminate points first
            if parm.fixedPoint
                AU = A(:,1:3*nPoints-1);
                AP = A(:,3*nPoints:end);
            else
                AU = A(:,1:3*nPoints);
                AP = A(:,(3*nPoints+1):end);
            end
            UU = AU'*AU + lambda*speye(size(AU,2));
            invUU = invertUU(UU,parm.fixedPoint);
            slask = AU'*AP;
            slask = invUU*slask;
            slask2 = AP'*AU;
            slask = slask2*slask;
            C = AP'*AP - slask;
            if size(C,2) == 0
                ATA_diag = 0;
            else
                ATA_diag = spdiags(sqrt(abs(diag(C))),0,sparse(size(C,2),size(C,2)));
            end
            C = (C+lambda*ATA_diag);
            slask = AU'*B;
            slask = invUU*slask;
            slask2 = AP'*AU;
            c = slask2*slask - AP'*B;
            if parm.verbosity>=2
                fprintf('\t\tSolving mod-Newton system.\n')
            end
            dR = C\c;
            dU = -invUU*(AU'*(AP*dR+B));
            d = [dU; dR];
        end
        [dpointvar, dcamvar] = restoreFixedVars(d, nPoints, nCams, ...
            parm.fixedCams, parm.fixedPoint, parm.fixAllPoints);
        [Pnew,Unew] = update_var(dpointvar, dcamvar, parm.P, parm.U);
        r = calcResidual(Pnew,Unew,parm);
        resnew = r'*r;
        if resnew > res
            lambda = lambda*2;
            if parm.verbosity>=2
                fprintf('\t%f\t%f', resnew, lambda)
            end
        end
    end
    if lambda > parm.minLambda
        lambda = max(parm.minLambda,lambda/1.5);
    end
    parm.U = Unew;
    parm.P = Pnew;
    if parm.verbosity>=2
        fprintf('%d\t%f\t%f',i,resnew,lambda)
    end
    if useProgressBar
        progressBar(100*i/iter);
    end
    %Check that residual shrinks "fast enough"
    if (res-resnew)/resnew < parm.minResChange
        if parm.verbosity>=2
            fprintf('\nResidual is not decreasing fast enough, exiting\n'); 
        end
        break
    end    
end

%Assign return parameters
P = parm.P;
U = parm.U;
for k = 1:length(P)
    P{k} = KK{k}*P{k};
end

  if parm.verbosity>=2
      fprintf('\nDone.\n')
  end
end

function tf = isCameraMatrixCellArray(arg)
    tf = iscell(arg) && all(2 == cellfun('ndims', arg)) && ...
        all(cellfun('isreal', arg)) && ...
        all(3==cellfun(@(x) size(x,1), arg)) && ...
        all(12==cellfun(@numel, arg));
end

function tf = isRelPoseArray(P,arg)
    N = numel(P);
    j = [arg.j];
    tf = isAbsPoseArray(P,arg) && isfield(arg, 'j') && all(isaninteger(j,1,N));
end

function tf = isPointTracks(P,arg)
    N = numel(P);
    tf = isstruct(arg) && isfield(arg, 'pointnr') && ...
        isfield(arg, 'points') && isfield(arg, 'index') && ...
        numel(arg.points) == N && numel(arg.index) == N && ...
        all(cellfun(@(x) all(isaninteger(x,1,arg.pointnr)), arg.index));
end

function tf = isAbsPoseArray(P,arg)
    N = numel(P);
    i = [arg.i];
    tf = isstruct(arg) && isfield(arg, 'i') && ...
        isfield(arg, 'M') && isfield(arg, 'w') && all(isaninteger(i,1,N));
end

function tf = isCamArray(P,arg)
    N = length(P);
    tf = all(isaninteger(arg,1,N));
end

function tf = isaninteger(x,low,high)
    tf = isnumeric(x) & x >= low & ...
        x <= high & round(x) == x;
end

function B = calcResidual(Pnew,Unew,parms)
    nPoints = parms.tracks.pointnr;
    B = [];
    if ~isempty(Unew)
        B = setupLinSystemCamera(Pnew,Unew,parms.tracks,...
            parms.camWeight,B);
    end
    if ~isempty(parms.absPos)
        B = setupLinSystemAbsPos(Pnew,parms.absPos,...
            nPoints,parms.posWeight,B);
    end
    if ~isempty(parms.relPoses)
        B = setupLinSystemRelCam(Pnew,parms.relPoses, ...
            nPoints, B);
    end
    if ~isempty(parms.absPoses)
        B = setupLinSystemAbsCam(Pnew,parms.absPoses, ...
            nPoints, B);
    end
end

function [B, A] = setupLinSystem(parms)
    P = parms.P;
    nCams = length(P);
    nPoints = parms.tracks.pointnr;
    nVars = 3*nPoints + 6*nCams;
    A = sparse(0,nVars);
    B = [];
    if ~isempty(parms.U)
        [B,A] = setupLinSystemCamera(P,parms.U,parms.tracks,...
            parms.camWeight,B,A);
    end
    if ~isempty(parms.absPos)
        [B,A] = setupLinSystemAbsPos(P,parms.absPos,...
            nPoints,parms.posWeight,B,A);
    end
    if ~isempty(parms.relPoses)
        [B,A] = setupLinSystemRelCam(P,parms.relPoses, ...
            nPoints, B, A);
    end
    if ~isempty(parms.absPoses)
        [B,A] = setupLinSystemAbsCam(P,parms.absPoses, ...
            nPoints, B, A);
    end
    
    %Remove columns corresponding to fixed variables
    %First the fixed cameras
    camindex = parms.fixedCams;
    if ~isempty(camindex)
        index = setdiff(1:length(P),camindex);
        keepcols = [(index-1)*6+1;(index-1)*6+2;(index-1)*6+3;(index-1)*6+4;(index-1)*6+5;(index-1)*6+6];
        keepcols = keepcols(:)';
        A = A(:,[1:3*nPoints (3*nPoints+keepcols)]);
    end
    %Fix points (first or all)
    if parms.fixAllPoints
        A = A(:,nPoints*3+1:end);        
    elseif parms.fixedPoint
        A = A(:,2:end);
    end
end

