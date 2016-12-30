function i=fitconic(points,normals);
% i=fitconic(points,normals);
% INPUT:
%   points - image point positions in homogeneous coordinates (3xn matrix).
%   normals - error normals for each point (3xn matrix)
% OUTPUT:
%   i - imagedata object
% If no normals are specified, the error is assumed to be isotropic.
%
% A future version i=fitconic(points,stddevs,normals);
% could incorporate known standard deviation in estimated point
% positions.
%   stddevs - estimated standard deviation for each image point (1xn matrix)

if nargin<=1,
  [un,L]=fitconic1(points);
else
  [un,L]=fitconic2(points,normals);
end
if isempty(un),
  i=[];
else
  i=imagedata;
  i=addconics(i,un,L);
end

function [un,L,n]=fitconic1(points);
% [un,L,n]=fitconic1(points);
% INPUT:
%   points - image point positions in homogeneous coordinates (3xn matrix).
% OUTPUT:
%   c - line parameters
%   L - cholesky factorisation of inverse of covariance of c.
%
% A future version [un,L,n]=fitconic1(points,stddevs);
% could incorporate known standard deviation in estimated point
% positions.
%   stddevs - estimated standard deviation for each image point (1xn matrix)

% Normalise points
  points=pflat(points);

  sc=1; %rescale points
  K=diag([1/sc 1/sc 1]);
  %K=eye(3);
  kpoints=K*points;
  % There is an error in the newton-raphson method.
  % Make a vector of all weights if you want to correct for
  % known standard deviation.
  % weights = [stddevs;stddevs]; weights = weights(:)/sc;

% Rough estimate
  M=[];
  dC=zeros(3,3,6);
  E=eye(6);
  for i=1:6, dC(:,:,i)=v2m(E(:,i)); end;
  for i=1:size(kpoints,2);
   M=[M; m2v2(kpoints(:,i)*kpoints(:,i)')'];
  end;
  [U,S,V]=svd(M);
  u=V(:,6);

  %Convert to ellipse parameters z, a, b, alpha
  A   = [u(1) u(2); u(2) u(3)];
  bb  = [2*u(4); 2*u(5)]; 
  c   = u(6);
  [Q D] = eig(A);
  det   = D(1,1)*D(2,2);
  if (det <= 0),
    warning('Degenerate ellipse');
    un=[];L=[];
  else
%    z = [0;0];
%    a = 1; b = 1; alpha = 0;
%  else 
    bs    = Q'*bb;
    alpha = atan2(Q(2,1), Q(1,1));
    zs    = -(2*D)\bs;  
    z     = Q*zs;
    h     = -bs'*zs/2-c;
    a     = sqrt(h/D(1,1));
    b     = sqrt(h/D(2,2));
%  end

  pv = [z;a;b;alpha];

 % Initialize phi (representation of points on conics closest to kpoints)
  X = kpoints(1:2,:)';
  c = cos(alpha); s = sin(alpha);
  Q = [c -s; s c];
  % compute initial approximations for phi_i
  du  = Q'*[ X(:,1)-z(1) X(:,2)-z(2)]';
  phi = angle(du(1,:)/a  + sqrt(-1)*du(2,:)/b)';

 % Newton-Raphson iteration
 for i=1:15,
  % Calculate residual Y
  s = sin(alpha);
  c = cos(alpha);
  Q = [c -s; s c];

  Xs = X*Q;
  zs = Q'*z;
  Y = [Xs(:,1)-zs(1)-a*cos(phi); Xs(:,2)-zs(2)-b*sin(phi)];
%  Y = Y ./ weights; % Correction for knonw standard deviation
%  norm(Y)
  res = Y;
  m = size(X,1);

  %% form Jacobian J
    S = sin(phi);  
    C = cos(phi);
    s = sin(alpha);
    c = cos(alpha);
    A = [-b*S C zeros(size(phi))  c*ones(size(phi)) s*ones(size(phi))];
    B = [ a*C zeros(size(phi)) S -s*ones(size(phi)) c*ones(size(phi))];

    [cg, sg] = rot_cossin (-a*S, b*C);
    G = sparse ([diag(cg), -diag(sg); diag(sg), diag(cg)]);
    Y = G*Y;

    D = diag (- a*S.*cg - b*C.*sg);
    J = [[D; zeros(m,m)], G*[A; B]];
%    J = [[D; zeros(m,m)], G*[A; B]] ./ (weights*ones(1,m+5));
%    Correction for knonw standard deviation.
%    OBS! This line is not correct yet!

    % Optional qr factorisation
    %[qq, J(m+1:2*m, m+1:m+5)] = qr(J(m+1:2*m, m+1:m+5));    
    %Y(m+1:m+5, :) = qq(:,1:5)'*Y(m+1:2*m,:);
    %h = J(1:m+5, 1:m+5)\Y(1:m+5);
    h = J\Y;

%   norm(h)
   phi = phi + h(1:m);
   alpha = alpha + h(m+1);
   a = a + h(m+2);
   b = b + h(m+3);
   z = z + h((m+4):(m+5));
  end;

  estimstd = sqrt((Y'*Y) / (size(Y,1)-size(Y,1)/2-5));
  Fx = J'*J;
  Cx = estimstd^2*inv(Fx);
  Cpv = Cx((m+1):(m+5),(m+1):(m+5));
  pv =[alpha;a;b;z];

  %Change to dual conic coordinates
  [D,u,dudpv]=ellipsparam2dualconic(pv);
  Cu = dudpv*Cpv*dudpv';

  % Correct for K
  up = m2v(inv(K)*v2m(u)*inv(K'));
  E=eye(6);
  dupdu=zeros(6,6);
  for i=1:6,
   dupdu(:,i) = -m2v( inv(K)'*v2m(E(:,i))*inv(K') );
  end;
  Cup = dupdu*Cu*dupdu';

  % PROJICERA NER up
  n=up/norm(up);
  un= up/(up'*n);
  dundup = (eye(size(up,1))/(up'*n) - up*n'/(up'*n)^2);
  Cun = dundup*Cup*dundup';

  % Cholesky factorisering
  [U,S,V]=svd(Cun);
  U=U(:,1:5);
  S=S(1:5,1:5);
  L=inv(sqrtm(S))*U';
  n=un;

end %not degenerate ellipse



function [c, s] = rot_cossin (x, y);
%ROT_COSSIN     Givens rotation angles
%
% [c, s] = rot_cossin (x, y);
% returns cos and sin vectors for Givens-rotation matrix
% which rotates y to zero.
%
% x, y: vectors
% c(i), s(i): [c(i) -s(i); s(i) c(i)]*[x(i); y(i)] == [..; 0]
  
  m = size(x,1);
  c = zeros(m,1); s = zeros(m,1);
  for i=1:m,
    if (abs(y(i)) > abs(x(i))),
      cot = -x(i)/y(i); si = 1/sqrt(1+cot^2); co = si*cot;
    else
      tan = -y(i)/x(i); co = 1/sqrt(1+tan^2); si = co*tan;
    end
    s(i) = si; c(i) = co;
  end 

%end % rot_cossin

function [C,u,dudpv]=ellipsparam2dualconic(paramvector);

alfa = paramvector(1);
a = paramvector(2);
b = paramvector(3);
z = paramvector(4:5);


D = diag([a^2 b^2 -1]);
T1 = [cos(alfa) -sin(alfa) 0; sin(alfa) cos(alfa) 0; 0 0 1];
T2 = [eye(2) z; 0 0 1];
C = T2*T1*D*T1'*T2';
dT2dz1 = [0 0 1;0 0 0;0 0 0];
dT2dz2 = [0 0 0;0 0 1;0 0 0];
dT1dalfa = [-sin(alfa) -cos(alfa) 0;cos(alfa) -sin(alfa) 0; 0 0 0];
dDda = diag([2*a 0 0]);
dDdb = diag([0 2*b 0]);

u=m2v(C);
dudpv = [ ...
m2v(T2*dT1dalfa*D*T1'*T2' + T2*T1*D*dT1dalfa'*T2') ...
m2v(T2*T1*dDda*T1'*T2') ...
m2v(T2*T1*dDdb*T1'*T2') ...
m2v(dT2dz1*T1*D*T1'*T2' + T2*T1*D*T1'*dT2dz1') ...
m2v(dT2dz2*T1*D*T1'*T2' + T2*T1*D*T1'*dT2dz2') ...
];


function v=m2v2(m);
% M2V v=m2v(m) returns vector of symmetric conic/quadric matrix

if size(m,1)==4,
   v=[m(1,1),2*m(1,2),m(2,2),2*m(1,3),2*m(2,3),m(3,3),...
	2*m(1,4),2*m(2,4),2*m(3,4),m(4,4)]';
elseif size(m,1)==3,
   v=[m(1,1),2*m(1,2),m(2,2),2*m(1,3),2*m(2,3),m(3,3)]';
else
  error('Wrong matrix-dimension');
end


function [un,L,n]=fitconic2(points,normals);
% [un,L,n]=fitline2(points,normals);
% INPUT:
%   points - image point positions in homogeneous coordinates (3xn matrix).
%   normals - error normals for each point (3xn matrix)
% OUTPUT:
%   c - line parameters
%   L - cholesky factorisation of inverse of covariance of c.
%
%
% A future version [un,L,n]=fitconic2(points,stddevs,normals);
% could incorporate known standard deviation in estimated point
% positions.
%   stddevs - estimated standard deviation for each image point (1xn matrix)

  sc=500; %rescale points

  K=diag([1/sc 1/sc 1]);
  kpoints=K*points;
  M=[];
  dC=zeros(3,3,6);
  E=eye(6);
  for i=1:6, dC(:,:,i)=v2m(E(:,i)); end;
  for i=1:size(kpoints,2);
   M=[M; m2v2(kpoints(:,i)*kpoints(:,i)')'];
  end;
  [U,S,V]=svd(M);
  u=V(:,6);
  ss=zeros(1,size(kpoints,2));
  for k=1:5;
   u=u/norm(u);
   C=v2m(u);
   for kk=1:5;
    idealpoints=kpoints+normals*diag(ss);
    %diag(idealpoints'*C*idealpoints)'
    deltass=-min(max((diag(idealpoints'*C*idealpoints))./ ...
          (diag(normals'*C*idealpoints)+...
           diag(idealpoints'*C*normals)),-0.02),0.02)'; % KOLLA DESSA GRÄNSER
    ss = ss+deltass;
    normss=norm(deltass);
   end;
   idealpoints=kpoints+normals*diag(ss);
   for i=1:6,
    dss1dC(:,i) = -diag(idealpoints'*dC(:,:,i)*idealpoints)./ ...
          (diag(normals'*C*idealpoints)+...
           diag(idealpoints'*C*normals)) + ...
              (diag(normals'*dC(:,:,i)*idealpoints)+...
           diag(idealpoints'*dC(:,:,i)*normals)).* ...
           diag(idealpoints'*C*idealpoints)./ ...
          (diag(normals'*C*idealpoints)+...
           diag(idealpoints'*C*normals)).^2;
    end;
%   u(1)=u(1)-litet;
%   u(2)=u(2)+litet;
%   ss1=ss';
%   [(ss1-ss0)/litet dss1dC(:,2)]
   [U,S,V]=svd(dss1dC);
   du = V(:,1:5)*inv(S(1:5,1:5))*U(:,1:5)'*(ss');
   u=u-du;
   normdu=norm(du);
   stddev=sqrt( (ss*ss')/(size(ss,2)-5) );
  end;
  %normss=norm(deltass)
  %normdu=norm(du)
  %stddev=sqrt( (ss*ss')/(size(ss,2)-5) );
                   % Estimate standard deviation in normal direction
                   % from residuals.
  [U,S,V]=svd(dss1dC);
  M= V(:,1:5)*inv(S(1:5,1:5))*U(:,1:5)';
  Cu = M*diag(stddev^2*ones(size(ss)))*M';
  % transformera till
  E=eye(6);
  C1 = K'*v2m(u)*K;
  C2 = inv( K'*v2m(u)*K );
  up = m2v(inv( K'*v2m(u)*K ));
  for i=1:6,
   dupdu(:,i) = -m2v( C2*K'*v2m(E(:,i))*K*C2 );
  end;
%  up0 = m2v(inv( K'*v2m(u)*K ));
%  up1 = m2v(inv( K'*v2m(u+litet*E(:,i))*K ));
%  [(up1-up0)/litet dupdu(:,i)];
  Cup = dupdu*Cu*dupdu';
  % PROJICERA NER up
  n=up/norm(up);
  un= up/(up'*n);
  dundup = (eye(size(up,1))/(up'*n) - up*n'/(up'*n)^2);
  Cun = dundup*Cup*dundup';
  [U,S,V]=svd(Cun);
  U=U(:,1:5);
  S=S(1:5,1:5);
  L=inv(sqrtm(S))*U';
  n=un;

% END OF fitconic2

