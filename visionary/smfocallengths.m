function [str mot] = smfocallengths(imseq, index)
% [str mot] = smfocallengths(imseq, index)
% Computes euclidean reconstruction from two views assuming two unknown
% focallengths and unity skew, aspect ratio & zero principal point.
% Input:
%   imseq - cell list of imagedata objects
%   index - reference to 2 views (optional, default = [1 2])
% Output:
%   mot - computed motion object
%   str - computed structure object
% NB: At least 8 points must be visible in both images

if(nargin < 2)
    index = [1 2];
end

if(length(index) ~= 2)
    error('smfocallengths requires exactly two views');
end

[str mot] = smproj8p(imseq, index);
F = ptobi(mot);

p1 = getpoints(imseq{1});
p2 = getpoints(imseq{2});
[f1, f2] = focallengths_inner(F, [0;0], [0;0]);

K1 = diag([f1 f1 1]);
K2 = diag([f2 f2 1]);

% compute essential matrix
E = K2 * F * K1;

% make sure that E is essential
[u s v] = svd(E);
s(end) = 0;
s([1 5]) = mean(s([1 5]));
%E = u*s*v';

% construct cameras
P = eye(3, 4);
w = [0 -1  0
    1  0  0
    0  0  1];

% factorize essential matrix to get cameras
% [u s v] = svd(E); (already done)
t = u(:, 3);
PP{1} = [u*w*v' t];
PP{2} = [u*w*v' -t];
PP{3} = [u*w'*v' t];
PP{4} = [u*w'*v' -t];

% triangulate a point to resolve factorization ambiguity.
for i = 1 : 4
    mm = motion({P, PP{i}});
    str = intsec2views(mm,...
        {imagedata([], p1(:,1)), imagedata([], p2(:,1))});
    X = getpoints(str);
    x = PP{i}*X;
    [x(end, :) X(3, :)]
    if(all([x(end, :) X(3, :)] > 0))
        % Depths are positive so this is the right one.
        P2 = PP{i};
        P1 = P;
        break;
    end
end

% uncalibrate cameras
P1 = K1*P1;
P2 = K2*P2;

mot = motion({P1, P2});

% triangulate
str = intsec2views(mot, imseq);

function [f1,f2]=focallengths_inner(F,pp1,pp2);
%focal lengths from F given pp1,pp2
%NB: Degenerate when principal axes intersect or equivalently
%      pp1 and pp2 are corresponding points <=> pp1'*F*pp2=0

pp1 = pextend(pp1);
pp2 = pextend(pp2);

%epipoles
e1=null(F');
e2=null(F);

E1=skew(e1);
E2=skew(e2);

I2=diag([1,1,0]);

numer1 = -pp2'*E2*I2*F'*pp1*pp1'*F*pp2;
denom1 = pp2'*E2*I2*F'*I2*F*pp2;
f1= sqrt(numer1/denom1);

numer2 = -pp1'*E1*I2*F*pp2*pp2'*F'*pp1;
denom2 = pp1'*E1*I2*F*I2*F'*pp1;
f2 = sqrt(numer2/denom2);