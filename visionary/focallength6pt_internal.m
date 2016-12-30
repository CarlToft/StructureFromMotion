function [sol, condition] = focallength6pt_internal(poly_eqs)
% for internal use with focallength6pt

% Remove redundant equations
[u s v] = svd(poly_eqs);
poly_eqs = v';
poly_eqs = poly_eqs(1:32,:);

% Eliminate
poly_eqs_elim = poly_eqs(:,1:32)\poly_eqs;

condition = cond(poly_eqs(:,1:32));
% Create action matrix
% Begin with the modulo matrix
mod_matrix = zeros(15,47);
mod_matrix(:,33:47) = eye(15);

% Complete the modulo matrix
mod_matrix(:,1:32) = -poly_eqs_elim(1:32,33:47)';

% Now construct the action matrix for l1
% These are the indices that are reached when the basis monomials are
% multiplied with l1.
ind_l1 = [20:24 28:31 33 34 38:40 44];
% Create the big action matrix
Ahat_l1 = zeros(47,15);
Ahat_l1(ind_l1,:) = eye(15);
% Reduce to the quotient space
A_l1 = mod_matrix*Ahat_l1;
% Calculate the eigen vectors of the transposed action matrix
[v1, d1] = eig(A_l1');
% Extract the solution for l1, l2 and p
sol = pflat(v1);
sol = sol(12:14, :);