function s=intsecconics(m,imseq,quadrics,option);
% INTSECCONICS s=intsecconics(m,imseq,quadrics,option) intersects quadrics from conics
%   with linear method
% INPUT:
%   m - motion object
%   imseq - cell array of imagedata objects
%   quadrics - (optional) specifies quadrics to be intersected. Otherwise all quadrics.
%            If 'quadrics' is a structure object, already existing quadrics are copied.
%   option:
%     'nocoordtrans' - No coordinate transformation
%     'conic' - Quadrics are diskquadrics, i.e. conics in a plane
%  
% OUTPUT:
%   s - structure object containing reconstructed 3D quadrics
% The results may be inaccurate.

if nargin<=2 | isempty(quadrics)
  quadrics = 1:size(imseq{1},3);
end

if nargin<=3
  option = [];
end
if isa(quadrics,'structure');
  nbrquadrics = size(quadrics,3);
  if nbrquadrics>0,
    U = getquadrics(quadrics);
  else
    nbrquadrics=size(imseq{1},3);
    U = NaN*ones(10,nbrquadrics);
  end
else
  nbrquadrics=length(quadrics);
  U = NaN*ones(10,nbrquadrics);
end

nbrimages=length(imseq);

for qq=1:nbrquadrics;

 if ~finite(U(1,qq)) %already existing?
  jj=0;
  M=[];
  for i=1:nbrimages;
   P=getcameras(m,i);
   c=getconics(imseq{i},qq);
   if finite(P(1)) & finite(c(1)),
    if strmatch('nocoordtrans',option),
     T=eye(3);
    else
     T=getnormtrans(imseq{i});
    end
    c=m2v(T*v2m(c)*T');
    conicP=coniccamera(T*P);
    z0=zeros(jj*6,1);
    z1=zeros(6,jj);
    M=[M,z0;conicP/norm(conicP,'fro'),z1,psphere(c)];
    jj=jj+1;
   end
  end

  if jj>2,
    [u,ss,v]=svd(M);
    if strmatch('conic',option);

% Old method. Find 'la' such that det(v1+la*v2)=0.
%     v1 = v(1:10,jj+10);
%     v2 = v(1:10,jj+9);
%     la=calclambdaconic(v1,v2);
%     L = v1 + la * v2;

      v1 = v(1:10,jj+10);
      [u,ss,v]=svd(v2m(v1));
      ss(4,4)=0;
      U(:,qq) = m2v(u*ss*v');
    else
      U(:,qq) = v(1:10,jj+10);
    end
  end
 end %already existing
end

s=structure([],[],U);


%end of main function

function la=calclambdaconic(L,K);
%L "best" quadric, K "second best"
% find best conic

K1=K(1);K2=K(2);K3=K(3);K4=K(4);K5=K(5);
K6=K(6);K7=K(7);K8=K(8);K9=K(9);K10=K(10);
L1=L(1);L2=L(2);L3=L(3);L4=L(4);L5=L(5);
L6=L(6);L7=L(7);L8=L(8);L9=L(9);L10=L(10);

%calc. polynom

      p4 = K7*K7*K5*K5-K1*K3*K9*K9-K4*K4*K3*K10-2.0*K4*K8*K7*K5-K7*K7*K3*K6+K1* ...
K3*K6*K10+2.0*K2*K8*K7*K6-2.0*K2*K8*K4*K9-K1*K5*K5*K10+K4*K4*K8*K8+K2*K2*K9*K9- ...
K2*K2*K6*K10+2.0*K1*K5*K8*K9+2.0*K2*K5*K4*K10+2.0*K4*K3*K7*K9-K1*K8*K8*K6-2.0* ...
K2*K5*K7*K9;

      tmp = 2.0*L2*K8*K7*K6-K1*L3*K9*K9+2.0*K4*K3*L7*K9+2.0*K4*L3*K7*K9-2.0*K4* ...
K8*K7*L5-2.0*K4*K8*L7*K5-2.0*K4*K3*L4*K10+2.0*K4*K3*K7*L9-K4*K4*K3*L10-2.0*L2* ...
K8*K4*K9-2.0*L2*K5*K7*K9+2.0*L2*K5*K4*K10-2.0*K2*L8*K4*K9+2.0*K2*L8*K7*K6+2.0* ...
K2*K8*K7*L6+2.0*K2*K8*L7*K6-2.0*K2*K8*K4*L9-2.0*K2*L5*K7*K9-2.0*K2*K8*L4*K9+2.0 ...
*K2*K5*K4*L10+2.0*K2*K5*L4*K10-2.0*K2*L2*K6*K10-2.0*K2*K5*K7*L9-2.0*K2*K5*L7*K9 ...
+2.0*K2*L5*K4*K10+2.0*L1*K5*K8*K9-L1*K5*K5*K10-L1*K8*K8*K6;
      p3 = tmp+2.0*K2*L2*K9*K9-K2*K2*K6*L10-K2*K2*L6*K10+2.0*K1*L5*K8*K9+2.0*K1* ...
K5*L8*K9+L1*K3*K6*K10+2.0*K2*K2*K9*L9+2.0*K4*K8*K8*L4-2.0*K1*K3*K9*L9+K1*L3*K6* ...
K10-2.0*K1*K5*L5*K10+2.0*K1*K5*K8*L9-2.0*K1*K8*L8*K6-K1*K5*K5*L10-K1*K8*K8*L6- ...
L1*K3*K9*K9+K1*K3*K6*L10+K1*K3*L6*K10-2.0*K7*K3*L7*K6-2.0*L4*K8*K7*K5-2.0*K4*L8 ...
*K7*K5+2.0*L4*K3*K7*K9-K7*K7*L3*K6+2.0*K7*K7*K5*L5-K7*K7*K3*L6-K4*K4*L3*K10+2.0 ...
*K4*K4*K8*L8+2.0*K7*K5*K5*L7;

      tmp = K4*K4*L8*L8+K2*K2*L9*L9+L2*L2*K9*K9+K8*K8*L4*L4-K1*K3*L9*L9-2.0*K4* ...
L3*L4*K10+2.0*K4*L3*K7*L9+2.0*K4*L3*L7*K9+4.0*K4*K8*L4*L8-2.0*K4*K8*L7*L5-2.0* ...
K4*K3*L4*L10+2.0*L2*K8*K7*L6+2.0*L2*K8*L7*K6-2.0*L2*L8*K4*K9+2.0*L2*L8*K7*K6 ...
-2.0*L2*K8*K4*L9-2.0*L2*L5*K7*K9-2.0*L2*K8*L4*K9+2.0*L2*L5*K4*K10+2.0*L2*K5*L4* ...
K10-2.0*L2*K5*K7*L9-2.0*L2*K5*L7*K9+2.0*K2*L8*K7*L6+2.0*K2*L8*L7*K6+2.0*L2*K5* ...
K4*L10+2.0*K2*K8*L7*L6-2.0*K2*L8*K4*L9-2.0*K2*L8*L4*K9-2.0*K2*K8*L4*L9-2.0*K2* ...
L5*L7*K9+2.0*K2*L5*L4*K10+K7*K7*L5*L5+2.0*K2*L5*K4*L10-2.0*K2*L5*K7*L9-2.0*K2* ...
L2*K6*L10-2.0*K2*L2*L6*K10+4.0*K2*L2*K9*L9+2.0*K2*K5*L4*L10-2.0*K2*K5*L7*L9+2.0 ...
*L1*L5*K8*K9;
      p2 = tmp-2.0*L1*K8*L8*K6-2.0*L1*K5*L5*K10+2.0*L1*K5*K8*L9+L1*K3*K6*L10+L1* ...
K3*L6*K10-2.0*L1*K3*K9*L9+L1*L3*K6*K10+K5*K5*L7*L7-L1*L3*K9*K9-L1*K5*K5*L10-L1* ...
K8*K8*L6+2.0*K1*K5*L8*L9+2.0*K1*L5*K8*L9+2.0*K1*L5*L8*K9-2.0*K1*K8*L8*L6+2.0*L1 ...
*K5*L8*K9-K2*K2*L6*L10-L2*L2*K6*K10+K1*L3*K6*L10+K1*L3*L6*K10-2.0*K1*L3*K9*L9 ...
-2.0*K1*K5*L5*L10-K1*L5*L5*K10-K1*L8*L8*K6+K1*K3*L6*L10-2.0*K7*K3*L7*L6-2.0*K7* ...
L3*L7*K6+4.0*K7*K5*L7*L5-2.0*L4*L8*K7*K5-2.0*L4*K8*K7*L5+2.0*L4*K3*K7*L9+2.0*L4 ...
*K3*L7*K9+2.0*L4*L3*K7*K9-2.0*K4*L8*K7*L5-2.0*K4*L8*L7*K5-2.0*L4*K8*L7*K5-K7*K7 ...
*L3*L6-K3*L7*L7*K6+2.0*K4*K3*L7*L9-K4*K4*L3*L10-K3*L4*L4*K10;

      tmp = 2.0*L2*L8*K7*L6+K1*L3*L6*L10+2.0*K1*L5*L8*L9-K1*L3*L9*L9-K1*L5*L5* ...
L10-K1*L8*L8*L6+2.0*L2*L8*L7*K6-2.0*K4*L3*L4*L10+2.0*K4*L3*L7*L9-2.0*K4*L8*L7* ...
L5+2.0*L2*K8*L7*L6-2.0*L2*L8*K4*L9-2.0*L2*L8*L4*K9-2.0*L2*L5*K7*L9-2.0*L2*L5*L7 ...
*K9-2.0*L2*K8*L4*L9-2.0*L2*K5*L7*L9+2.0*L2*L5*K4*L10+2.0*L2*L5*L4*K10+2.0*L2*K5 ...
*L4*L10+2.0*K2*L5*L4*L10+2.0*L1*K5*L8*L9+2.0*L1*L5*K8*L9+2.0*L1*L5*L8*K9-2.0*L1 ...
*K8*L8*L6+2.0*K2*L8*L7*L6-2.0*L1*K5*L5*L10+L1*K3*L6*L10;
      p1 = tmp+L1*L3*K6*L10+L1*L3*L6*K10-L1*L5*L5*K10-L1*L8*L8*K6+2.0*K2*L2* ...
L9*L9-2.0*L1*L3*K9*L9-2.0*K2*L2*L6*L10-2.0*K2*L5*L7*L9-2.0*K2*L8*L4*L9-L2*L2*K6 ...
*L10-L2*L2*L6*K10+2.0*L2*L2*K9*L9+2.0*K4*L4*L8*L8-L1*K3*L9*L9-2.0*L4*K8*L7*L5 ...
-2.0*L4*L8*K7*L5-2.0*L4*L8*L7*K5+2.0*L4*L3*L7*K9+2.0*L4*K3*L7*L9+2.0*L4*L3*K7* ...
L9-2.0*K7*L3*L7*L6-K3*L7*L7*L6-L3*L7*L7*K6+2.0*K5*L7*L7*L5-K3*L4*L4*L10-L3* ...
L4*L4*K10+2.0*K8*L4*L4*L8+2.0*K7*L7*L5*L5;

      p0 = L1*L3*L6*L10-L1*L3*L9*L9-L1*L5*L5*L10+2.0*L1*L5*L8*L9-L1*L8*L8*L6- ...
L2*L2*L6*L10+L2*L2*L9*L9+2.0*L2*L5*L4*L10-2.0*L2*L5*L7*L9-2.0*L2*L8*L4*L9+2.0* ...
L2*L8*L7*L6-L3*L4*L4*L10+2.0*L4*L3*L7*L9+L4*L4*L8*L8-2.0*L4*L8*L7*L5-L3*L7*L7* ...
L6+L7*L7*L5*L5;

r=roots([p4,p3,p2,p1,p0]);

rindex=find(real(r)==r);
%rindex=find(abs(imag(r))<1e-4);

if length(rindex)>0,
 rabs=abs(real(r(rindex)));
 rmin=find(min(rabs)==rabs);
 la = real(r(rindex(rmin(1))));
else
 la = 0; %not possible
end

