function Evec = calibrated_fivepoint(Q1, Q2, init)
% Evec = calibrated_fivepoint(Q1, Q2, init)
%
% If init == 0 it will not check if there is a compiled version of 
% calibrated_fivepoint_helper.c
%
% Q1 end Q2 are correspoding image points.
%
% Code to verify that it works: 
% Q1 = rand(3,5);
% Q2 = rand(3,5);
% Evec = calibrated_fivepoint( Q1,Q2);
% for i=1:size(Evec,2)
%   E = reshape(Evec(:,i),3,3);
%   % Check determinant constraint! 
%   det( E)
%   % Check trace constraint
%   2 *E*transpose(E)*E -trace( E*transpose(E))*E
%   % Check reprojection errors
%   diag( Q1'*E*Q2)
% end
%
% PS: Please note that due to varying standards of which is Q1 and Q2
% it is very possible that you get essential matrices which are 
% the transpose of what your expected. 

if nargin < 3
    init = 1;
end

if init
    try
        calibrated_fivepoint_helper();
    catch le;
        if strcmp(le.identifier,'MATLAB:UndefinedFunction')
            mex calibrated_fivepoint_helper.c
        end
    end
end

Evec = calibrated_fivepoint_gb(Q1,Q2);
