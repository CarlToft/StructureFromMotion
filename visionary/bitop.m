function m=bitop(F);
% function m=bitop(F);
% Determines two possible camera matrices from the fundamental matrix F
%
% Input: F=fundamental matrix
% Output: [P1,P2]=two camera matrices normalised such that
%          P1=[I 0] P2(1,:)=[0 0 0 1]
%         The 2 cameras are returned as MOTION m 

[E12,E21]=bitoone(F);
E=E12./E12(1);
P1=[1 0 0 0;0 1 0 0;0 0 1 0];
P2=[0 0 0 1;F(1,3) F(2,3) F(3,3) E(2);-F(1,2) -F(2,2) -F(3,2) E(3)];

m=motion({P1,P2});
