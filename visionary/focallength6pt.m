function mot = focallength6pt(imseq)
% mot = focallength(imseq) solves the 6pt-2view minimal case for cameras
% with a common unknown focallength. Has 15 solutions in general but only
% real solutions with positive depths for all points are returned. This
% means typically 1-4 solutions.
%
% Input:
%   imseq - cellist of imagedata. 2 views with 6 points each.
%
% Output:
%   mot - cellist of motion objects containing two cameras each.

% get first 6 points in each image
ims_old = imseq;
x1 = getpoints(imseq{1});
x1 = x1(:, 1:6);
imseq{1} = imagedata([], x1);
x2 = getpoints(imseq{2});
x2 = x2(:, 1:6);
imseq{2} = imagedata([], x2);

% normalize data for numerical performance
scale = max(abs([x1(:); x2(:)]));
H = diag([1/scale 1/scale 1]);
imseq{1} = changecsystem(imseq{1}, H);
imseq{2} = changecsystem(imseq{2}, H);

% generate equations
[eqs, F_parts] = focallength6pt_create_eqs(imseq);

% solve equations
% [sol, condition] = focallength6pt_svd_internal(eqs);
% [sol, condition] = classic_camera_varf(eqs);
[sol, condition] = focallength6pt_internal(eqs);
    
% reformat output
for k = 1 : size(sol, 2)
    F{k} = F_parts{1} + sol(1,k)*F_parts{2} + sol(2,k)*F_parts{3};
    f(k) = 1 / sqrt(sol(3,k));
end

% pick only real solutions
F = F(f==real(f));
f = f(f==real(f));

% construct cameras
w = [0 -1  0
     1  0  0
     0  0  1];
kk = 0;
for k = 1 : length(f)
    if(~isreal(F{k}))
        continue;
    end
    kk = kk + 1;
    % factorize essential matrix to get cameras
    K = diag([f(k) f(k) 1]);
    E = K'*F{k}'*K; % we transpose F to comply with the H-Z definition of F.
    [u s v] = svd(E);
    t = u(:, 3);
    PP{1} = K*[u*w*v' t];
    PP{2} = K*[u*w*v' -t];
    PP{3} = K*[u*w'*v' t];
    PP{4} = K*[u*w'*v' -t];
    P = K*eye(3, 4);
    
    % triangulate a point to resolve factorization ambiguity.
    best = -inf;
    for i = 1 : 4
        X = getpoints(intsec2views(motion({P, PP{i}}), imseq));
        x = PP{i}*X;
        X3 = X(3,:);
        X3 = X3 / norm(X3);
        x3 = x(3,:);
        x3 = x3 / norm(x3);
        score = sum(x3) + sum(X3);
        if(score > best)
            best = score;
            bestcam = i;
        end
    end

    P2 = H\PP{bestcam};
    P1 = H\P;
    mot{kk} = motion({P1, P2});
end