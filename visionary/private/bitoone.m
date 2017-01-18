function [E12,E21]=bitoone(F);

% function [E12,E21]=bitoone(F);
% Determines the one forms between image 1 and 2 from
% the fundamental matrix between images 1 and 2,
% with a correct scale between them.
%
% Input  F=fundamental matrix between images 1 and 2
% Output [E12,E21]=the one forms between image 1 and 2 and
%                  betweem image 2 and 1 respectively

AF=[...
det(F([2,3],[2,3])) det(F([2,3],[3,1])) det(F([2,3],[1,2]));...
det(F([3,1],[2,3])) det(F([3,1],[3,1])) det(F([3,1],[1,2]));...
det(F([1,2],[2,3])) det(F([1,2],[3,1])) det(F([1,2],[1,2]))];

[U,S,V]=svd(AF);

E21=sqrt(S(1,1)).*U(:,1)';
E12=sqrt(S(1,1)).*V(:,1)';




