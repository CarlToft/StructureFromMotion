function s=imconic(points,normals,stddevs);
% IMCONIC/IMCONIC constructor
% s=imconic(points,normals,stddevs);
%

if nargin == 0,
  s.points = [];
  s.normals = [];
  s.stddevs = [];
  s.idealpoints = [];
  s.u = [];
  s.L = [];
  s.n = [];
  s = class(s,'imconic');
elseif isa(points,'imconic');
  ps = points.points;
  ns = points.normals;
  stds = points.stddevs;
  s=imconic(ps,ns,stds);
else
%  keyboard;
  K=diag([1/500 1/500 1]);
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
           diag(idealpoints'*C*normals)),-0.02),0.02)';
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
   stddev=std(ss);
  end;
  normss=norm(deltass)
  normdu=norm(du)
  stddev=std(ss);
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

  s.points  = points;
  s.normals = normals;
  s.stddevs = stddevs;
  s.idealpoints = inv(K)*idealpoints;
  s.u = un;
  s.L = L;
  s.n = n;
  s = class(s,'imconic');
end;
