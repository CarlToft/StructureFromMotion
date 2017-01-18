function [T,Kout]=flexnonlinear(mot,pp,option);
% function [T,K]=flexnonlinear(mot,pp,option);
% Computes projective to Euclidean upgrade nonlinearly
% optimizing min(sum(norm(Ki*Ki')-norm(Pi*Om*Pi'))).
% Pollefey's method for initialization.
% Input:
%   mot - projective motion
%   pp - approximate principal point (optional)
%   option:
%     'autocalib=XXXXX' - Autocalibration
%       specify for each of focal length, aspect ratio, skew,
%       principal point x, principal point y,
%       if the parameter is
%        1 known and nominal,
%        2 unknown but constant in the sequence,
%        3 unknown and varying.
%     'nolinearinit' - No linear initialization (Tinit=eye(4))
% Output:
%   T - 4 x 4 coordinate transformation for Euclidean upgrade
%       Use changecsystem(mot,T) for upgrade!
%   Kout - cell list of calibration matrices
% Default: pp=[0,0]'; and constant intrinsic parameters

if nargin<2,
  pp=[0,0]';
end;
linearinit=1;
caliboptions=2*ones(1,5);

if nargin>=3,
  if strmatch('nolinearinit',option);
    linearinit=0;
  end
    
  if strmatch('autocalib=',option),
    q=strmatch('autocalib=',option);
    strautocalib = option{q}(11:length(option{q}));
    caliboptions=[str2num(strautocalib(1)) str2num(strautocalib(2)) ...
                  str2num(strautocalib(3)) str2num(strautocalib(4)) ...
		  str2num(strautocalib(5))];
  end
end

if linearinit,
  T0=flexlinear(mot,pp);
  T0=T0/norm(T0);
  iT0=inv(T0);
else
  T0=eye(4);
  iT0=eye(4);
end

nbrcameras=size(mot);
Afull=iT0(1:3,1:3);
[A,R]=rq(Afull);
c=R*iT0(4,1:3)'/A(3,3);
A=A/A(3,3);

%Om=[A*A', A*c;c'*A',c'*c]
%Om=inv(T0)*diag([1,1,1,0])*inv(T0)';
%Om=Om/norm(Om);

Kmean=zeros(3,3);
P=zeros(3,4*nbrcameras);
Kout=cell(1,nbrcameras);

for i=1:nbrcameras;
  Pi=getcameras(mot,i);
  [K,R]=rq(Pi*iT0);
  K=K/K(3,3);K([1,5])=abs(K([1,5]));
  Kmean=Kmean+K;
end
Kmean=Kmean/nbrcameras;
keyboard;

index=[1,4,5,7,8];

params=[A(index)';c(:)];
Kvnom=[1,1,0,0,0];
Kvmean=[Kmean(1,1),Kmean(2,2)/Kmean(1,1),Kmean(1,2)/Kmean(1,1),...
		Kmean(1,3),Kmean(2,3)];

focalmean=Kvmean(1);

for i=1:nbrcameras,
  Pi=getcameras(mot,i);
  [K,R]=rq(Pi*iT0);
  K=K/K(3,3);
  Kv=[K(1,1),K(2,2)/K(1,1),K(1,2)/K(1,1),K(1,3),K(2,3)];

%  for j=1:5,
%    if caliboptions(j)==2 & i==1,
%      params=[params;Kvmean(j)];
%    elseif caliboptions(j)==3,
%      params=[params;Kv(j)];
%    end
%  end
%  P(:,4*i-3:4*i)=Pi/norm(Pi);
%end


  for j=1:5,
    if caliboptions(j)==1,
      Kv(j)=Kvnom(j);
    elseif caliboptions(j)==2,
      Kv(j)=Kvmean(j);
      if i==1,
        params=[params;Kv(j)/focalmean];
      end
    elseif caliboptions(j)==3,
      params=[params;Kv(j)/focalmean];
    end
  end;
  K0=[Kv(1),Kv(1)*Kv(3),Kv(4);0,Kv(1)*Kv(2),Kv(5);0,0,1];
  Pi=diag([1/focalmean,1/focalmean,1])*Pi;
  P(:,4*i-3:4*i)=Pi/norm(Pi);
end


params=absconstopt(P,params,caliboptions);

A=eye(3);A(index)=params(1:5);
c=params(6:8);
T=inv([A,zeros(3,1);c',1]);


indexp=9;
Kvmean=zeros(1,5);

for i=1:nbrcameras,
  Kv=[1,1,0,0,0];
  for j=1:5,
    if caliboptions(j)==2,
      if i==1,
        Kvmean(j)=params(indexp);
	indexp=indexp+1;
      end;
      Kv(j)=Kvmean(j);
    elseif caliboptions(j)==3,
      Kv(j)=params(indexp);
      indexp=indexp+1;
    end
  end
  Kout{i}=diag([focalmean,focalmean,1])*[Kv(1),Kv(1)*Kv(3),Kv(4);0,Kv(1)*Kv(2),Kv(5);0,0,1];
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,f0] = absconstopt(P,A,caliboptions);


[f0, Idiff, Idiffgrad] = absconstres(P,A,caliboptions);
dA=inf;
counter=0;
while log(norm(dA)+eps)/log(10)>-15 & counter<1000;
%  log(norm(dA)+eps)/log(10)
  S = Idiffgrad'*Idiffgrad;
  dA=-inv(S+0.0000001*eye(size(S,1)))*Idiffgrad'*Idiff;
%  Anew=A+reshape(dA,size(A))';
  Anew=A+dA;
  [f1, Idiff, Idiffgrad] = absconstres(P,Anew,caliboptions);
  while f1>f0 & counter<1000,
%    f0
%    f1
%    disp('Det har tar vi en gang till');
    dA=dA/2;
%    Anew=A+reshape(dA,size(A))';
    Anew=A+dA;
    [f1, Idiff, Idiffgrad] = absconstres(P,Anew,caliboptions);
    counter=counter+1;
  end;

  if f1<=f0,
    f0=f1;
    A = Anew;
  end;
  counter=counter+1;
end;

function [f,Idiff,Idiffgrad]=absconstres(P,params,caliboptions);

nbrcameras=size(P,2)/4;
Idiff=zeros(6*nbrcameras,1);
nbrvar=8+sum(caliboptions==2)+sum(caliboptions==3)*nbrcameras;
Idiffgrad=zeros(6*nbrcameras,nbrvar);


index=[1,4,5,7,8];
index2=[index,9];
A=eye(3);A(index)=params(1:5);
c=params(6:8);




%new
Om=[A*A',A*c;c'*A',c'*c];


dOm1=zeros(4,4);dOm1(1,1)=2*A(1,1);dOm1(1,4)=c(1);dOm1(4,1)=c(1);
dOm2=zeros(4,4);dOm2(1,1)=2*A(1,2);dOm2(1,2)=A(2,2);dOm2(2,1)=A(2,2);
dOm2(1,4)=c(2);dOm2(4,1)=c(2);
dOm3=zeros(4,4);dOm3(1,2)=A(1,2);dOm3(2,1)=A(1,2);dOm3(2,2)=2*A(2,2);
dOm3(2,4)=c(2);dOm3(4,2)=c(2);
dOm4=zeros(4,4);dOm4(1,1)=2*A(1,3);dOm4(1,2)=A(2,3);dOm4(2,1)=A(2,3);
dOm4(1,3)=A(3,3);dOm4(3,1)=A(3,3);dOm4(1,4)=c(3);dOm4(4,1)=c(3);
dOm5=zeros(4,4);dOm5(1,2)=A(1,3);dOm5(2,1)=A(1,3);dOm5(2,2)=2*A(2,3);
dOm5(2,3)=A(3,3);dOm5(3,2)=A(3,3);dOm5(2,4)=c(3);dOm5(4,2)=c(3);
%dOm6=zeros(4,4);dOm6(1,3)=A(1,3);dOm6(3,1)=A(1,3);dOm6(2,3)=A(2,3);
%dOm6(3,2)=A(2,3);dOm6(3,3)=2*A(3,3);dOm6(3,4)=c(3);dOm6(4,3)=c(3);
dOm7=zeros(4,4);dOm7(1,4)=A(1,1);dOm7(4,1)=A(1,1);dOm7(4,4)=2*c(1);
dOm8=zeros(4,4);dOm8(1,4)=A(1,2);dOm8(4,1)=A(1,2);dOm8(2,4)=A(2,2);dOm8(4,2)=A(2,2);dOm8(4,4)=2*c(2);
dOm9=zeros(4,4);dOm9(1,4)=A(1,3);dOm9(4,1)=A(1,3);dOm9(2,4)=A(2,3);dOm9(4,2)=A(2,3);dOm9(3,4)=A(3,3);dOm9(4,3)=A(3,3);dOm9(4,4)=2*c(3);

indexp=9;
Kvmean=zeros(1,5);
Kvmeanindex=zeros(1,5);

for i=1:nbrcameras,

%first term

  Kv=[1,1,0,0,0];
  Kvindex=zeros(1,5);

  for j=1:5,
    if caliboptions(j)==2,
      if i==1,
	Kvmeanindex(j)=indexp;
        Kvmean(j)=params(indexp);
	indexp=indexp+1;
      end
      Kv(j)=Kvmean(j);
      Kvindex(j)=Kvmeanindex(j);
    elseif caliboptions(j)==3,
      Kv(j)=params(indexp);
      Kvindex(j)=indexp;
      indexp=indexp+1;
    end
  end
  K=[Kv(1),Kv(1)*Kv(3),Kv(4);0,Kv(1)*Kv(2),Kv(5);0,0,1];
  KK=K*K';
  N1=sum(KK(:).^2);

  for j=1:5,
    if Kvindex(j)>0,
      dKK=zeros(3,3);
      switch j,
        case 1
	  dKK(1,1)=2*Kv(1)*(1+Kv(3)^2);
	  dKK(1,2)=2*Kv(1)*Kv(2)*Kv(3);
	  dKK(2,1)=2*Kv(1)*Kv(2)*Kv(3);
	  dKK(2,2)=2*Kv(1)*Kv(2)^2;
        case 2
	  dKK(1,2)=Kv(1)^2*Kv(3);
	  dKK(2,1)=Kv(1)^2*Kv(3);
	  dKK(2,2)=2*Kv(1)^2*Kv(2);
	case 3
	  dKK(1,1)=2*Kv(1)^2*Kv(3);
	  dKK(1,2)=Kv(1)^2*Kv(2);
	  dKK(2,1)=Kv(1)^2*Kv(2);
	case 4
	  dKK(1,1)=2*Kv(4);
	  dKK(1,2)=Kv(5);
	  dKK(2,1)=Kv(5);
	  dKK(1,3)=1;
	  dKK(3,1)=1;
        case 5
	  dKK(1,2)=Kv(4);
	  dKK(2,1)=Kv(4);
	  dKK(2,2)=2*Kv(5);
	  dKK(2,3)=1;
	  dKK(3,2)=1;
      end
      df=dKK/sqrt(N1)-KK/N1^(3/2)*sum(sum(KK.*dKK));
      Idiffgrad(6*i-5:6*i,Kvindex(j))=df(index2)';
    end
  end



%second term
  p=P(:,4*i-3:4*i);
  pop=p*Om*p';
  N2=sum(pop(:).^2);

  res=KK/sqrt(N1)-pop/sqrt(N2);
  Idiff(6*i-5:6*i)=res(index2);

  pdop1=p*dOm1*p';
  dfO1=pdop1/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop1));
  Idiffgrad(6*i-5:6*i,1)=-dfO1(index2)';

  pdop2=p*dOm2*p';
  dfO2=pdop2/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop2));
  Idiffgrad(6*i-5:6*i,2)=-dfO2(index2)';

  pdop3=p*dOm3*p';
  dfO3=pdop3/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop3));
  Idiffgrad(6*i-5:6*i,3)=-dfO3(index2)';

  pdop4=p*dOm4*p';
  dfO4=pdop4/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop4));
  Idiffgrad(6*i-5:6*i,4)=-dfO4(index2)';

  pdop5=p*dOm5*p';
  dfO5=pdop5/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop5));
  Idiffgrad(6*i-5:6*i,5)=-dfO5(index2)';

%  pdop6=p*dOm6*p';
%  dfO6=pdop6/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop6));
%  Idiffgrad(6*i-5:6*i,6)=-dfO6(index2)';

  pdop7=p*dOm7*p';
  dfO7=pdop7/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop7));
  Idiffgrad(6*i-5:6*i,6)=-dfO7(index2)';

  pdop8=p*dOm8*p';
  dfO8=pdop8/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop8));
  Idiffgrad(6*i-5:6*i,7)=-dfO8(index2)';

  pdop9=p*dOm9*p';
  dfO9=pdop9/sqrt(N2)-pop/N2^(3/2)*sum(sum(pop.*pdop9));
  Idiffgrad(6*i-5:6*i,8)=-dfO9(index2)';
end;


f=Idiff'*Idiff;


%if 0,
%
% slask=0.001;
% param0=params;
% Idiffgr0=Idiffgrad;
% Idiff0=Idiff;
% params=param0;
% params(9)=params(9)+slask;
%
% [Idiffgr0(:,1)  (Idiff-Idiff0)/slask]
%
%end;
