
function [t,R,R2]=invessentialm(E,u1,u2);
% [t,R,R2]=invessentialm(E,u1,u2) - decompose the essential matrix
% into a unit direction t and a rotation matrix so that
% E=essentialm(t,R) (up to scale)
% INPUT:
%   E - essential matrix
%   u1,u2 - 3xn homogenous image points (optional)
%
% OUTPUT
%   t,R - translation rotation
%   R2 - second choice of rotation. Not needed if u1,u2 are specified 
%        since this solution can be discarded by chirality
P=[0  1 0;-1 0 0; 0 0 -1]; % perm. matrix
P2=[0  1 0;-1 0 0; 0 0 1]; % perm. matrix
Q=[0 -1 0; 1 0 0; 0 0 0]; % cross matrix

[uE,sE,vE]=svd(E);

s1=sE(1,1);
s2=sE(2,2);
if abs (s1-s2)>1e-3,
  error('Not essential matrix');
else
  % s(3,3) should be zero
  s3 = s2; % to make S full rang
  % F = u*P*Q*inv(u*P)*u*P*diag([s1 s2 s3])*v'

  T= uE*P*Q*inv(uE*P);
  t=[ T(3,2) T(1,3) T(2,1)]';
  R=uE*P*diag([1 1 det(uE)*det(vE)])*vE';
  R=sign(det(R))*R;
  R2=uE*P2*diag([1 1 det(uE)*det(vE)])*vE';
  R2=sign(det(R2))*R2;
  
  if nargin>1,
        % Choose the one that gives as much positive depth vectors as
        % possible.
        nr_of_points=size(u1,2);
        P1=[eye(3) zeros(3,1)];
        P2=[R' -R'*t];
        % Calculate the structure U.
        for i=1:nr_of_points,
            M=[P1 u1(:,i) zeros(3,1); P2 zeros(3,1) u2(:,i)];
            % One ought to check that the intersection is well defined.
            [u,s,v]=svd(M);
            %diag(s)'
            U(1:4,i)=v(1:4,6);
            %l1(1,i)=-v(5,6);
            %l2(1,i)=-v(6,6);
        end;
        [U,aa]=pflat(U);
        [u1t,l1]=psphere(P1*U);
        [u2t,l2]=psphere(P2*U);
        o1=sign(sum(u1t.*u1));
        o2=sign(sum(u2t.*u2));
        % the vectors o1 and o2 indicate whether the points in 
        % image 1 and 2 have positive depth or not. Ideally
        % all elements of both o1 and o2 should be positive.
        if sum(o1.*o2)<0,
            R=R2; %If they have different signs change R from R1 to R2.
            P2=[R' -R'*t];
            % Calculate the structure U.
            for i=1:nr_of_points,
                M=[P1 u1(:,i) zeros(3,1); P2 zeros(3,1) u2(:,i)];
                % One ought to check that the intersection is well defined.
                [u,s,v]=svd(M);
                %diag(s)'
                U(1:4,i)=v(1:4,6);
                %l1(1,i)=-v(5,6);
                %l2(1,i)=-v(6,6);
            end;
            U=pflat(U);
            u1t=psphere(P1*U);
            u2t=psphere(P2*U);
            o1=sign(sum(u1t.*u1));
            o2=sign(sum(u2t.*u2));
        end;
        if sum(o2)<0,
            t=-t; % If they have negative sign change t from t to -t.
        end;
        % Now we have done our best at getting positive depth in the camera
        % equation. Calculate once more the camera matrices P1 and P2,
        % and the reprojected points u1 and u2 and the sign of the depths o1
        % and o2.
  end
end

