function [U,P,s] = krot4_feas_LP_sparse(u,A,tol,min_depth,max_depth)

numcams = length(A);

% remove points not visible in two views.
vis = zeros(1,u.pointnr);
for i=1:length(u.points)
    imvis = zeros(1,u.pointnr);
    imvis(u.index{i}) = 1;
    vis = vis + imvis;
end
vis0 = vis;
vis = vis >= 2;
for i = 1:length(u.points)
    pp = NaN*ones(3,u.pointnr);
    pp(:,u.index{i}) = u.points{i};
    pp = pp(:,vis);
    u.index{i} = find(isfinite(pp(1,:)));
    u.points{i} = pp(:,isfinite(pp(1,:)));
end
u.pointnr = sum(vis);

%Set up the matrices for the problem
[a,a0,b,b0,c,c0] = gen_krot(u,A);
%Solve
[Linfsol,s] = LinfSolverfeas(a,a0,b,b0,c,c0,tol,min_depth,max_depth);
%Reshape the solution to the correct form
[U,P] = form_str_mot(u,A,Linfsol);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,P] = form_str_mot(u,A,sol)
numpts = u.pointnr;
numcams = length(A);

U = reshape(sol(1:(3*numpts)), [3 numpts]);
U = pextend(U);

tpart = sol(3*numpts+1:end);
P = cell(size(A));
P{1} = [A{1} [0 0 0]'];
for i=2:length(A)
    P{i} = [A{i} [tpart([(i-2)*3+1:(i-1)*3])]];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y,s] = LinfSolverfeas(a,a0,b,b0,c,c0,tol,min_depth,max_depth)
	
A1 = [-a-tol*c; a-tol*c; -b-tol*c; b-tol*c];
B1 = [a0+tol*c0; -a0+tol*c0; b0+tol*c0; -b0+tol*c0];

%slack
A1 = [A1 [-speye(size(a,1)); speye(size(a,1)); sparse(2*size(a,1),size(a,1))]];
A1 = [A1 [sparse(2*size(a,1),size(a,1)); -speye(size(a,1)); speye(size(a,1))]];
%L1 term (|s_i| < t_i minimize sum t_i)
A1 = [A1 sparse(size(A1,1),2*size(a,1))];
A1 = [A1; ...
    sparse(2*size(a,1),size(a,2)) -speye(2*size(a,1)) -speye(2*size(a,1));...
    sparse(2*size(a,1),size(a,2)) speye(2*size(a,1)) -speye(2*size(a,1))];
B1 = [B1; sparse(4*size(a,1),1)];

% Depth limitations
A2 = [-c sparse(size(a,1),4*size(a,1)); c sparse(size(a,1),4*size(a,1))];
B2 = [c0-min_depth; max_depth-c0];

A = [A1; A2];
buc = [B1; B2];

C = [sparse(size(a,2),1); sparse(2*size(a,1),1); ones(2*size(a,1),1)];

par = msklpopt(C,A,[],buc,[],[],[],'param');
par.param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; %Turn off simplex solver
par.param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 3000; %increase max interations
par.param.MSK_IPAR_LOG = 0; %increase max interations
[m,n] = size(a);
clear A1 A2 B1 B2 a a0 b b0 c c0;
result = msklpopt(C,A,[],buc,[],[],par.param,'minimize');

if result.rcode~=0,
%    error('ERROR IN MOSEK');
    disp('ERROR IN MOSEK');
end
Y = result.sol.itr.xx;
X = result.sol.itr.suc;

s = Y(n+2*m+1:end);
s = s(1:m)+s(m+1:2*m);
Y = Y(1:n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,a0,b,b0,c,c0] = gen_krot(u,A)
numvar = 3*u.pointnr+3*length(A);
numpts = u.pointnr;
numcams = length(A);
%Set up the problem
a = [];
b = [];
c = [];
for i = 1:numcams;
    R = A{i};
    p = NaN*ones(3,u.pointnr);
    p(:,u.index{i}) = u.points{i};
    
    visible_points = isfinite(p(1,:));
    numres = sum(visible_points);
    
    %First term
    %compute coefficients in from tor every term.
    ptind = find(visible_points');
    pointcoeff = p(1,visible_points)'*R(3,:)-ones(numres,1)*R(1,:);
    %row and column in the A-matrisen
    pointcol = [(ptind-1)*3+1 (ptind-1)*3+2 ptind*3];
    pointrow = [1:numres]'*[1 1 1];
    
    %Coeff. in front of the translation part of the camera
    tcoeff = [-ones(numres,1) zeros(numres,1) p(1,visible_points)'];
    tcol = ones(numres,1)*[numpts*3+[(i-1)*3+1:i*3]];
    trow = pointrow;
    
    %Data for A matrix
    data = [pointcoeff(:); tcoeff(:)];
    row = [pointrow(:);  trow(:)];
    col = [pointcol(:); tcol(:)];
    newa = sparse(row,col,data,numres,numvar);

    
    %Second term
    %coefficients for the 3D points.
    ptind = find(visible_points');
    pointcoeff = p(2,visible_points)'*R(3,:)-ones(numres,1)*R(2,:);
    %row and column of the B-matrix
    pointcol = [(ptind-1)*3+1 (ptind-1)*3+2 ptind*3];
    pointrow = [1:numres]'*[1 1 1];
    
    %Coeff for the translation
    tcoeff = [zeros(numres,1) -ones(numres,1) p(2,visible_points)'];
    tcol = ones(numres,1)*[numpts*3+[(i-1)*3+1:i*3]];
    trow = pointrow;
    
    %Data for the B-matrix.
    data = [pointcoeff(:); tcoeff(:)];
    row = [pointrow(:);  trow(:)];
    col = [pointcol(:); tcol(:)];
    newb = sparse(row,col,data,numres,numvar);

    %The demnominator
    %Coefficients for the 3D-points
    ptind = find(visible_points');
    pointcoeff = ones(numres,1)*R(3,:);
    
    pointcol = [(ptind-1)*3+1 (ptind-1)*3+2 ptind*3];
    pointrow = [1:numres]'*[1 1 1];
    
    %Coefficients for the translation
    tcoeff = [zeros(numres,1) zeros(numres,1) ones(numres,1)];
    tcol = ones(numres,1)*[numpts*3+[(i-1)*3+1:i*3]];
    trow = pointrow;
    
    %Data
    data = [pointcoeff(:); tcoeff(:)];
    row = [pointrow(:);  trow(:)];
    col = [pointcol(:); tcol(:)];
    newc = sparse(row,col,data,numres,numvar);

    if 0;
        UU = rand(3,numpts);
        UUU = UU(:);
        t = rand(3, numcams);
        tt = t(:);
        var = [UUU; tt];
        slask = R*UU + repmat(t(:,i),[1 size(UU,2)]);
        figure(1);plot((p(1,visible_points).*slask(3,visible_points)-slask(1,visible_points))'-newa*var)
        figure(2);plot((p(2,visible_points).*slask(3,visible_points)-slask(2,visible_points))'-newb*var)
        figure(3);plot(slask(3,visible_points)'-newc*var)
    end

    %Add to previous data.
    a = [a;newa];
    b = [b;newb];
    c = [c;newc];
end

%Select coordinate system such that 
%first camcenter = (0,0,0), scale is determined by depth conditions.

a = a(:,[1:numpts*3 (numpts*3+4):end]);
b = b(:,[1:numpts*3 (numpts*3+4):end]);
c = c(:,[1:numpts*3 (numpts*3+4):end]);
a0 = zeros(size(a,1),1);
b0 = zeros(size(b,1),1);
c0 = zeros(size(c,1),1);

if 0;
    R = A{i};
    UU = rand(3,numpts);
    UU(:,1) = 0;
    UUU = UU(:);
    t = rand(3, numcams);
    t(3,1) = 1;
    tt = t(:);
    var = [UUU; tt];
    slask = R*UU + repmat(t(:,i),[1 size(UU,2)]);
    figure(1);plot((p(1,visible_points).*slask(3,visible_points)-slask(1,visible_points))'-newa*var)
    figure(2);plot((p(2,visible_points).*slask(3,visible_points)-slask(2,visible_points))'-newb*var)
    figure(3);plot(slask(3,visible_points)'-newc*var)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);

function y = pextend(x)
y = [x; ones(1,size(x,2))];
