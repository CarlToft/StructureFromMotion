function [U,P,lambda] = modbundle_sparse(U,P,u,iter,lambda,goodindex,ddfactor)
% optimering av \sum_j ||f_ij(P_i,U_j)-u_j||_2^2 
% genom linearisering av f_ij(P_i,U_j)
% iter - antal iterationer
% lambda - dämpfaktor

if nargin<6,
    goodindex=[];
end
if nargin<7,
    ddfactor = 10;
end
schur = 1;

nbr=u.pointnr;
u.pointnr = u.pointnr + ddfactor*length(goodindex);
for jj=1:length(goodindex);
    U(:,end+[1:ddfactor])=U(:,goodindex(jj))*ones(1,ddfactor);
end

for i = 1:length(u.points);
    [KKK,R] = rq(P{i}(:,1:3));
    KK{i} = KKK;
    u.points{i} = pflat(inv(KKK)*u.points{i});
    [aa,bb]=intersect(u.index{i},goodindex);
    if length(bb)>0,
        for jj=1:length(bb);
            ind = find(goodindex==aa(jj));
            u.index{i}(end+[1:ddfactor])=nbr+ddfactor*(ind-1)+[1:ddfactor];
            u.points{i}(:,end+[1:ddfactor])=u.points{i}(:,bb(jj))*ones(1,ddfactor);
        end
    end
    P{i} = inv(KKK)*P{i};
end

fprintf('Iter:\t Error:\t lambda:\n');
res = compute_res(P,U,u);                 
fprintf('%d\t%f\t%f',0,res,lambda);
for i = 1:iter
    [A,B] = setup_lin_system(P,U,u);  
    res = compute_res(P,U,u);
    if ~schur
        C = (A'*A+lambda*speye(size(A,2),size(A,2)));
        c = A'*B;
        fprintf('\t\t\t\t\t\t\tSolving mod-Newton system.');
        d = -C\c;
    else
        %%%% Eliminera punkter %%%%%%%%%%
        AU = A(:,1:3*u.pointnr-1);
        AP = A(:,3*u.pointnr:end);
        UU = AU'*AU + lambda*speye(size(AU,2));
        invUU = invertUU(UU);
        slask = AU'*AP;
        slask = invUU*slask;
        slask2 = AP'*AU;
        slask = slask2*slask;
        C = lambda*speye(size(slask)) + AP'*AP - slask;
        slask = AU'*B;
        slask = invUU*slask;
        slask2 = AP'*AU;
        c = slask2*slask - AP'*B;
        fprintf('\t\t\t\t\t\t\tSolving mod-Newton system.');
        dR = C\c;
        dU = -invUU*(AU'*(AP*dR+B));
        d = [dU; dR];
    end    
    fprintf('\tDone.\n');
    [Pnew,Unew] = update_var(d,P,U);
    resnew = compute_res(Pnew,Unew,u);
    while resnew > res
        lambda = lambda*2;
        if ~schur
            C = (A'*A+lambda*speye(size(A,2),size(A,2)));
            c = A'*B;
            fprintf('\t\t\t\t\t\t\tSolving mod-Newton system.');
            d = -C\c;
        else
            %%%% Eliminera punkter %%%%%%%%%%
            AU = A(:,1:3*u.pointnr-1);
            AP = A(:,3*u.pointnr:end);
            UU = AU'*AU + lambda*speye(size(AU,2));
            invUU = invertUU(UU);
            slask = AU'*AP;
            slask = invUU*slask;
            slask2 = AP'*AU;
            slask = slask2*slask;
            C = lambda*speye(size(slask)) + AP'*AP - slask;
            slask = AU'*B;
            slask = invUU*slask;
            slask2 = AP'*AU;
            c = slask2*slask - AP'*B;
            fprintf('\t\t\t\t\t\t\tSolving mod-Newton system.');
            dR = C\c;
            dU = -invUU*(AU'*(AP*dR+B));
            d = [dU; dR];
        end
        [Pnew,Unew] = update_var(d,P,U);
        resnew = compute_res(Pnew,Unew,u);
    end
    if lambda > 0.001
        lambda = lambda/1.5;
    end
    U = Unew;
    P = Pnew;
    fprintf('%d\t%f\t%f',i,resnew,lambda);
end
fprintf('\n');

for i = 1:length(P);
    P{i} = KK{i}*P{i};
end
U=U(:,1:nbr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = compute_res(P,U,u)
res = 0;
for i = 1:length(P);
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = pflat(u.points{i});
    vis = isfinite(uu(1,:));
    res = res + ...
        sum(((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - uu(1,vis)).^2) + ...
        sum(((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - uu(2,vis)).^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pnew,Unew] = update_var(d,P,U)
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];

dpointvar = [0; d(1:(3*size(U,2)-1))];
dpointvar = reshape(dpointvar, size(U(1:3,:)));
dcamvar = [0;0;0;0;0;0;d(3*size(U,2):end)];
dcamvar = reshape(dcamvar,[6 length(P)]);

Unew = pextend(U(1:3,:) + dpointvar);

Pnew = cell(size(P));
for i=1:length(P);
    R0 = P{i}(:,1:3);
    t0 = P{i}(:,4);
    R = expm(Ba*dcamvar(1,i) + Bb*dcamvar(2,i) + Bc*dcamvar(3,i))*R0;
    t = t0 + dcamvar(4:6,i);
    Pnew{i} = [R t];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = setup_lin_system(P,U,u)
numpts = size(U,2);
%Bas för tangentplanet till rotationsmångfalden.
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];

da1 = cell(size(P));
db1 = cell(size(P));
dc1 = cell(size(P));
dt11 = cell(size(P));
dt21 = cell(size(P));
dt31 = cell(size(P));
dU11 = cell(size(P));
dU21 = cell(size(P));
dU31 = cell(size(P));
da2 = cell(size(P));
db2 = cell(size(P));
dc2 = cell(size(P));
dt12 = cell(size(P));
dt22 = cell(size(P));
dt32 = cell(size(P));
dU12 = cell(size(P));
dU22 = cell(size(P));
dU32 = cell(size(P));


U = pflat(U);
U = U(1:3,:);

for i=1:length(P);
    %beräkna derivator för båda residualeran i alla bilder
    %a,b,c - rotations parametrar för kamera i
    %U1,U2,U3 - 3d punkt parametrar
    %t1,t2,t3 - translations parametrar kameran.
    vis = u.index{i};
    R0 = P{i}(:,1:3);
    %   t0 = repmat(P{i}(:,4),[1 size(u{i},2)]);
    t0 = P{i}(:,4);
    da1{i} =  (Ba(1,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(1,:)*U(:,vis)+t0(1))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Ba(3,:)*R0*U(:,vis));

    da2{i} = (Ba(2,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(2,:)*U(:,vis)+t0(2))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Ba(3,:)*R0*U(:,vis));
    
    db1{i} = (Bb(1,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(1,:)*U(:,vis)+t0(1))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Bb(3,:)*R0*U(:,vis));
    
    db2{i} = (Bb(2,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(2,:)*U(:,vis)+t0(2))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Bb(3,:)*R0*U(:,vis));
    
    dc1{i} = (Bc(1,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(1,:)*U(:,vis)+t0(1))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Bc(3,:)*R0*U(:,vis));

    dc2{i} = (Bc(2,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(2,:)*U(:,vis)+t0(2))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Bc(3,:)*R0*U(:,vis));
    
    dU11{i} = R0(1,1)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(1,:)*U(:,vis) + t0(1))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,1);

    dU12{i} = R0(2,1)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(2,:)*U(:,vis) + t0(2))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,1);
    
    dU21{i} = R0(1,2)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(1,:)*U(:,vis) + t0(1))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,2);

    dU22{i} = R0(2,2)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(2,:)*U(:,vis) + t0(2))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,2);
    
    dU31{i} = R0(1,3)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(1,:)*U(:,vis) + t0(1))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,3);
    
    dU32{i} = R0(2,3)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(2,:)*U(:,vis) + t0(2))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,3);

    dt11{i} = 1./(R0(3,:)*U(:,vis)+t0(3));
    dt12{i} = zeros(size(dt11{i}));
    
    dt21{i} = zeros(size(dt11{i}));
    dt22{i} = 1./(R0(3,:)*U(:,vis)+t0(3));
    
    dt31{i} = -(R0(1,:)*U(:,vis)+t0(1))./((R0(3,:)*U(:,vis)+t0(3)).^2);
    dt32{i} = -(R0(2,:)*U(:,vis)+t0(2))./((R0(3,:)*U(:,vis)+t0(3)).^2);
        
    if 0
        aa = 2*rand(1,1)-1;
        bb = 2*rand(1,1)-1;
        cc = 2*rand(1,1)-1;
        UU1 = 2*rand(1,1)-1;
        UU2 = 2*rand(1,1)-1;
        UU3 = 2*rand(1,1)-1;
        tt1 = 2*rand(1,1)-1;
        tt2 = 2*rand(1,1)-1;
        tt3 = 2*rand(1,1)-1;       

        f = [];
        vis = find(isfinite(u{i}(1,:)));
        kk = ceil(rand*length(vis));
        kk = vis(kk);
        f01 = (-u{i}(1,kk)+(R0(1,:)*U(:,kk)+ t0(1))./(R0(3,:)*U(:,kk)+t0(3)));
        f02 = (-u{i}(2,kk)+(R0(2,:)*U(:,kk)+t0(2))./(R0(3,:)*U(:,kk)+ t0(3)));
        
        f2 = [];
        for s = -0.1:0.001:0.1
            R = expm(s*(Ba*aa+Bb*bb+Bc*cc))*R0;
            t = t0+s*[tt1;tt2;tt3];
            UU = U+s*repmat([UU1;UU2;UU3], [1 size(U,2)]);
            f = [f; (u{i}(1,kk)-(R(1,:)*UU(:,kk)+t(1))./(R(3,:)*UU(:,kk)+ t(3)))^2 + ...
                    (u{i}(2,kk)-(R(2,:)*UU(:,kk)+t(2))./(R(3,:)*UU(:,kk)+ t(3)))^2];
            f2 = [f2; (f01 + s*(da1{i}(kk)*aa + db1{i}(kk)*bb + dc1{i}(kk)*cc + ...
                    dU11{i}(kk)*UU1 + dU21{i}(kk)*UU2 + dU31{i}(kk)*UU3 + ...
                    dt11{i}(kk)*tt1 + dt21{i}(kk)*tt2 + dt31{i}(kk)*tt3)).^2 + ...
                    (f02 + s*(da2{i}(kk)*aa + db2{i}(kk)*bb + dc2{i}(kk)*cc + ...
                    dU12{i}(kk)*UU1 + dU22{i}(kk)*UU2 + dU32{i}(kk)*UU3 + ...
                    dt12{i}(kk)*tt1 + dt22{i}(kk)*tt2 + dt32{i}(kk)*tt3)).^2];
        end
        %figure(1);plot(-0.1:0.001:0.1,[f f2(end:-1:1)]);
        figure(1);plot(-0.1:0.001:0.1,[f f2]);
    end
end
row = [];
col = [];
data = [];
resnum = 0;
B = [];
%fprintf('\tSetting up system, camera: ');
fprintf('\tSetting up system. ');
%Allokera hela vektorn först två residualer för varje bild punkt
resnr = 0;
for i = 1:length(P);
    resnr = resnr + size(u.points{i},2);
end
row = zeros(2*resnr*(6+3),1);
col = zeros(2*resnr*(6+3),1);
data = zeros(2*resnr*(6+3),1);
B = zeros(2*resnr,1);
lastentry = 0;
lastentryB = 0;
for i = 1:length(P);
    %fprintf('%d ',i);
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = pflat(u.points{i});
    vis = find(isfinite(uu(1,:)));
    
    %Första residualen:
    %3D-punkt parametrar:
    %U1-koeff
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [(vis-1)*3+1]';
    data(lastentry+[1:length(vis)]) = dU11{i}';
    lastentry = lastentry+length(vis);
    
    %U2-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; [(vis-1)*3+2]'];
    %data = [data; dU21{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [(vis-1)*3+2]';
    data(lastentry+[1:length(vis)]) = dU21{i}';
    lastentry = lastentry+length(vis);
    
    %U3-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; [vis*3]'];
    %data = [data; dU31{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [vis*3]';
    data(lastentry+[1:length(vis)]) = dU31{i}';
    lastentry = lastentry+length(vis);
    
    %Kameraparametrar
    %a-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+1)*ones(length(vis),1)];
    %data = [data; da1{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+1)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = da1{i}';
    lastentry = lastentry+length(vis);
    
    %b-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+2)*ones(length(vis),1)];
    %data = [data; db1{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+2)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = db1{i}';
    lastentry = lastentry+length(vis);
    
    %c-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+3)*ones(length(vis),1)];
    %data = [data; dc1{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+3)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dc1{i}';
    lastentry = lastentry+length(vis);
    
    %t_1-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+4)*ones(length(vis),1)];
    %data = [data; dt11{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+4)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt11{i}';
    lastentry = lastentry+length(vis);
    
    %t_2-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+5)*ones(length(vis),1)];
    %data = [data; dt21{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+5)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt21{i}';
    lastentry = lastentry+length(vis);
    
    %t_3-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+i*6)*ones(length(vis),1)];
    %data = [data; dt31{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+i*6)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt31{i}';
    lastentry = lastentry+length(vis);
    
    %Andra residualen:
    %3D-punkt parametrar:
    %U1-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; [(vis-1)*3+1]'];
    %data = [data; dU12{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [(vis-1)*3+1]';
    data(lastentry+[1:length(vis)]) = dU12{i}';
    lastentry = lastentry+length(vis);
    
    
    %U2-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; [(vis-1)*3+2]'];
    %data = [data; dU22{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [(vis-1)*3+2]';
    data(lastentry+[1:length(vis)]) = dU22{i}';
    lastentry = lastentry+length(vis);
    
    %U3-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; [vis*3]'];
    %data = [data; dU32{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [vis*3]';
    data(lastentry+[1:length(vis)]) = dU32{i}';
    lastentry = lastentry+length(vis);

    %Kameraparametrar
    %a-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+1)*ones(length(vis),1)];
    %data = [data; da2{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+1)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = da2{i}';
    lastentry = lastentry+length(vis);
    
    %b-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+2)*ones(length(vis),1)];
    %data = [data; db2{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+2)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = db2{i}';
    lastentry = lastentry+length(vis);
    
    %c-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+3)*ones(length(vis),1)];
    %data = [data; dc2{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+3)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dc2{i}';
    lastentry = lastentry+length(vis);
    
    %t_1-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+4)*ones(length(vis),1)];
    %data = [data; dt12{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+4)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt12{i}';
    lastentry = lastentry+length(vis);
    
    %t_2-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+5)*ones(length(vis),1)];
    %data = [data; dt22{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+5)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt22{i}';
    lastentry = lastentry+length(vis);
        
    %t_3-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+i*6)*ones(length(vis),1)];
    %data = [data; dt32{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+i*6)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt32{i}';
    lastentry = lastentry+length(vis);
        
    resnum = resnum+2*length(vis);
    
    %Konstant-termer
    btmp = zeros(2*length(vis),1);
    %Första residualen
    btmp(1:2:end) = (P{i}(1,:)*pextend(U(:,vis)))./(P{i}(3,:)*pextend(U(:,vis)))-uu(1,vis);
    %Andra residualen
    btmp(2:2:end) = (P{i}(2,:)*pextend(U(:,vis)))./(P{i}(3,:)*pextend(U(:,vis)))-uu(2,vis);
    %B = [B; btmp];
    B(lastentryB+[1:length(btmp)]) = btmp;
    lastentryB = lastentryB+length(btmp);
    
    if 0
        A = sparse(row,col,data);
        UU = 2*rand(3,numpts)-1;
        pvar = 2*rand(6,length(P))-1;
        var = [UU(:); pvar(:)];
        PP = pvar(:,i);
        vis = isfinite(u{i}(1,:));
        aa = pvar(1,i);
        bb = pvar(2,i);
        cc = pvar(3,i);
        tt1 = pvar(4,i);
        tt2 = pvar(5,i);
        tt3 = pvar(6,i);
        
        UU1 = UU(1,:);
        UU2 = UU(2,:);
        UU3 = UU(3,:);
        f1 = da1{i}*aa + db1{i}*bb + dc1{i}*cc + ...
            dU11{i}.*UU1 + dU21{i}.*UU2 + dU31{i}.*UU3 + ...
            dt11{i}*tt1 + dt21{i}*tt2 + dt31{i}*tt3;
        f2 = da2{i}*aa + db2{i}*bb + dc2{i}*cc + ...
            dU12{i}.*UU1 + dU22{i}.*UU2 + dU32{i}.*UU3 + ...
            dt12{i}*tt1 + dt22{i}*tt2 + dt32{i}*tt3;
        f12 = (A*var(1:size(A,2)));
        
        figure(1);plot(f1(vis));
        figure(2);plot(f12(1:2:end));
        figure(1);plot(f1(vis)-f12(1:2:end)');
        figure(2);plot(f2(vis)-f12(2:2:end)');
        
        A = sparse(row,col,data);
        dU = 2*rand(size(U(1:3,:)))-1;
        dP = 2*rand(6,length(P));
        var = [dU(:); dP(:)];
        f = [];
        f2 = [];
        slask = -0.01:0.0001:0.01;
        for lambda = slask;
            Unew = pextend(U(1:3,:) + lambda*dU);
            R0 = P{i}(:,1:3);
            t0 = P{i}(:,4);
            R = expm(lambda*(Ba*dP(1,i) + Bb*dP(2,i) + Bc*dP(3,i)))*R0;
            t = t0+lambda*dP(4:6,i);
            Pnew = [R t];
            vis = isfinite(u{i}(1,:));
            f = [f;sum(sum((pflat(Pnew*Unew(:,vis))-u{i}(:,vis)).^2))];
            f2 = [f2; sum((A*var(1:size(A,2))*lambda+B).^2)];
        end
        figure(1);plot(slask,[f f2]);
        
    end
end

A = sparse(row,col,data);
%lås koordinatsystemet
%Första kameran konstant
%och första koordinaten i första punkten konstant
A = A(:,[1:3*numpts 3*numpts+7:end]);
A = A(:,[2:end]);
fprintf('Done\n');







if 0
    for i = 1:length(P);
    %fprintf('%d ',i);
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = pflat(u.points{i});
    vis = find(isfinite(uu(1,:)));
    
    %Första residualen:
    %3D-punkt parametrar:
    %U1-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; [(vis-1)*3+1]'];
    data = [data; dU11{i}'];
    %U2-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; [(vis-1)*3+2]'];
    data = [data; dU21{i}'];
    %U3-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; [vis*3]'];
    data = [data; dU31{i}'];
    %Kameraparametrar
    %a-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+1)*ones(length(vis),1)];
    data = [data; da1{i}'];
    %b-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+2)*ones(length(vis),1)];
    data = [data; db1{i}'];
    %c-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+3)*ones(length(vis),1)];
    data = [data; dc1{i}'];
    %t_1-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+4)*ones(length(vis),1)];
    data = [data; dt11{i}'];
    %t_2-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+5)*ones(length(vis),1)];
    data = [data; dt21{i}'];
    %t_3-koeff
    row = [row; resnum+[1:2:2*length(vis)]'];
    col = [col; (3*numpts+i*6)*ones(length(vis),1)];
    data = [data; dt31{i}'];
    
    %Andra residualen:
    %3D-punkt parametrar:
    %U1-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; [(vis-1)*3+1]'];
    data = [data; dU12{i}'];
    %U2-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; [(vis-1)*3+2]'];
    data = [data; dU22{i}'];
    %U3-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; [vis*3]'];
    data = [data; dU32{i}'];
    %Kameraparametrar
    %a-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+1)*ones(length(vis),1)];
    data = [data; da2{i}'];
    %b-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+2)*ones(length(vis),1)];
    data = [data; db2{i}'];
    %c-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+3)*ones(length(vis),1)];
    data = [data; dc2{i}'];
    %t_1-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+4)*ones(length(vis),1)];
    data = [data; dt12{i}'];
    %t_2-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+(i-1)*6+5)*ones(length(vis),1)];
    data = [data; dt22{i}'];
    %t_3-koeff
    row = [row; resnum+[2:2:2*length(vis)]'];
    col = [col; (3*numpts+i*6)*ones(length(vis),1)];
    data = [data; dt32{i}'];
    resnum = resnum+2*length(vis);
    
    %Konstant-termer
    btmp = zeros(2*length(vis),1);
    %Första residualen
    btmp(1:2:end) = (P{i}(1,:)*pextend(U(:,vis)))./(P{i}(3,:)*pextend(U(:,vis)))-uu(1,vis);
    %Andra residualen
    btmp(2:2:end) = (P{i}(2,:)*pextend(U(:,vis)))./(P{i}(3,:)*pextend(U(:,vis)))-uu(2,vis);
    B = [B; btmp];
    
    if 0
        A = sparse(row,col,data);
        UU = 2*rand(3,numpts)-1;
        pvar = 2*rand(6,length(P))-1;
        var = [UU(:); pvar(:)];
        PP = pvar(:,i);
        vis = isfinite(u{i}(1,:));
        aa = pvar(1,i);
        bb = pvar(2,i);
        cc = pvar(3,i);
        tt1 = pvar(4,i);
        tt2 = pvar(5,i);
        tt3 = pvar(6,i);
        
        UU1 = UU(1,:);
        UU2 = UU(2,:);
        UU3 = UU(3,:);
        f1 = da1{i}*aa + db1{i}*bb + dc1{i}*cc + ...
            dU11{i}.*UU1 + dU21{i}.*UU2 + dU31{i}.*UU3 + ...
            dt11{i}*tt1 + dt21{i}*tt2 + dt31{i}*tt3;
        f2 = da2{i}*aa + db2{i}*bb + dc2{i}*cc + ...
            dU12{i}.*UU1 + dU22{i}.*UU2 + dU32{i}.*UU3 + ...
            dt12{i}*tt1 + dt22{i}*tt2 + dt32{i}*tt3;
        f12 = (A*var(1:size(A,2)));
        
        figure(1);plot(f1(vis));
        figure(2);plot(f12(1:2:end));
        figure(1);plot(f1(vis)-f12(1:2:end)');
        figure(2);plot(f2(vis)-f12(2:2:end)');
        
        A = sparse(row,col,data);
        dU = 2*rand(size(U(1:3,:)))-1;
        dP = 2*rand(6,length(P));
        var = [dU(:); dP(:)];
        f = [];
        f2 = [];
        slask = -0.01:0.0001:0.01;
        for lambda = slask;
            Unew = pextend(U(1:3,:) + lambda*dU);
            R0 = P{i}(:,1:3);
            t0 = P{i}(:,4);
            R = expm(lambda*(Ba*dP(1,i) + Bb*dP(2,i) + Bc*dP(3,i)))*R0;
            t = t0+lambda*dP(4:6,i);
            Pnew = [R t];
            vis = isfinite(u{i}(1,:));
            f = [f;sum(sum((pflat(Pnew*Unew(:,vis))-u{i}(:,vis)).^2))];
            f2 = [f2; sum((A*var(1:size(A,2))*lambda+B).^2)];
        end
        figure(1);plot(slask,[f f2]);
        
    end
end

end

function invUU = invertUU(UU)
% fixa så att alla block är 3 block
[i,j,data] = find(UU);
UU = sparse(i+1,j+1,data);
UU(1,1) = 1;

%plocka ut elementen
[i,j,data] = find(UU);
%much faster if I use full matrices here
A = full(sparse(i,mod(j-1,3)+1,data));
%A = sparse(i,mod(j-1,3)+1,data);
invUU = speye(size(UU,1));
[i,j,data] = find(invUU);
invA = full(sparse(i,mod(j-1,3)+1,data));
%invA = sparse(i,mod(j-1,3)+1,data);

%Gaussa
%Nedåt eliminering
for i = 1:3;
    %div:a med diagonal element
    invA(i:3:end,:) = invA(i:3:end,:)./repmat(A(i:3:end,i),[1 3]);
    A(i:3:end,:) = A(i:3:end,:)./repmat(A(i:3:end,i),[1 3]);
    for j = i+1:3
        %dra bort element
        invA(j:3:end,:) = invA(j:3:end,:) - repmat(A(j:3:end,i),[1 3]).*invA(i:3:end,:);
        A(j:3:end,:) = A(j:3:end,:) - repmat(A(j:3:end,i),[1 3]).*A(i:3:end,:);
    end
end
%Uppåt eliminering
for i = 3:-1:2
    for j = i-1:-1:1;
        invA(j:3:end,:) = invA(j:3:end,:) - repmat(A(j:3:end,i),[1 3]).*invA(i:3:end,:);
        A(j:3:end,:) = A(j:3:end,:) - repmat(A(j:3:end,i),[1 3]).*A(i:3:end,:);
    end
end
[i,j,data] = find(invA);
invUU = sparse(i(2:end)-1,j(2:end)+(ceil(i(2:end)/3)-1)*3-1,data(2:end));

function [r,q]=rq(a)
% RQ [r,q]=rq(a) factorises a such a=rq where r upper tri. and q unit matrix
% If a is not square, then q is equal q=[q1 q2] where q1 is unit matrix

[m,n]=size(a);
e=eye(m);
p=e(:,[m:-1:1]);
[q0,r0]=qr(p*a(:,1:m)'*p);

r=p*r0'*p;
q=p*q0'*p;

fix=diag(sign(diag(r)));
r=r*fix;
q=fix*q;

if n>m
  q=[q, inv(r)*a(:,m+1:n)];
end

function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);

function y = pextend(x)
y = [x; ones(1,size(x,2))];

