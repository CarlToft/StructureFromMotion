function q=rot2quat(R);
% function q=rot2quat(R);
% Rotation matrix R to unit quaternion q

[slask,index]=max(diag(R));
if index==1,
    ii=[1,2,3];
elseif index==2,
    ii=[2,3,1];
else
    ii=[3,1,2];
end
r=sqrt(1+sum(diag(R)));
q=zeros(4,1);
q(1)=r/2;
q(ii(1)+1)=(R(ii(3),ii(2))-R(ii(2),ii(3)))/2/r;
q(ii(2)+1)=(R(ii(1),ii(3))-R(ii(3),ii(1)))/2/r;
q(ii(3)+1)=(R(ii(2),ii(1))-R(ii(1),ii(2)))/2/r;
