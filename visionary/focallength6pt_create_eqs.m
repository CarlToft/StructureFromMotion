function [poly_eqs, F_part] = focallength6pt_create_eqs(imseq)
% this is an internal function used by focallength6pt.m

x1 = getpoints(imseq{1});
x2 = getpoints(imseq{2});

%
% Construct the Fundamental matrix
%

% Setup all constraints on the fundamental matrix
M = [x1(1,:).*x2(1,:);
    x1(1,:).*x2(2,:);
    x1(1,:).*x2(3,:);
    x1(2,:).*x2(1,:);
    x1(2,:).*x2(2,:);
    x1(2,:).*x2(3,:);
    x1(3,:).*x2(1,:);
    x1(3,:).*x2(2,:);
    x1(3,:).*x2(3,:);
    ]';

% Find the nullspace of M
[U,S,V] = svd(M,0);
FF = V(:,7:9);

% The three parts of the fundamental matrix
F_part{1} = reshape(FF(:,1),3,3)';
F_part{2} = reshape(FF(:,2),3,3)';
F_part{3} = reshape(FF(:,3),3,3)';

%
% Construct the equations
%

% Setup the monomials
l1 = multipol(1,[1 0 0]');
l2 = multipol(1,[0 1 0]');
p = multipol(1,[0 0 1]');

% Create the fundamental matrix
F = F_part{1}+l1*F_part{2}+l2*F_part{3};

% Create all equations according to figure 3 of stewes paper;
tmp = [F(1,1) F(1,2) F(1,3)*p;...
       F(2,1) F(2,2) F(2,3)*p;...
       F(3,1) F(3,2) F(3,3)*p] *...
      [F(1,1) F(2,1) F(3,1)*p;...
       F(1,2) F(2,2) F(3,2)*p;...
       F(1,3) F(2,3) F(3,3)*p];
   
eqs_tmp = 2*tmp*F-trace(tmp)*F; 

eqs = {};
for ii = 1:9
    eqs{ii} = eqs_tmp(ii);
end

% det(F) = 0
eqs{10} = det(F);
eqs{11} = eqs{10}*p;
eqs{12} = eqs{11}*p;

% Multiply all equations with p and p^2 to get enough equations
eqs_p = {};
eqs_p2 = {};
for ii = 1:length(eqs);
    eqs_p{ii} = eqs{ii}*p;
    eqs_p2{ii} = eqs_p{ii}*p;
end
eqs = [eqs eqs_p eqs_p2];

% Build the coefficent matrix
[poly_eqs mons] = polynomials2matrix(eqs);

% Change the order of the monomials to simplify things
poly_eqs(:,35:37) = poly_eqs(:,[37 35 36]);

% Remove the columns that will not be correctly eliminated
poly_eqs(:,[4 20 30]) = [];
