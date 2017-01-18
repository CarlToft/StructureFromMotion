function [dPdx,ppout,focalout,arout] = accalcdPdx(K,R,t,caliboptions,lockcsystem,pp,focal,ar);

focalout=[];
focalcnt=0;
arout=[];
arcnt=0;
ppout=[];
ppcnt=0;

b1=[0 1 0;-1 0 0;0 0 0];
b2=[0 0 -1;0 0 0;1 0 0];
b3=[0 0 0;0 0 1;0 -1 0];

nbrcam=size(K,1);
index=0;
dkindex=zeros(nbrcam,5);
for jj=1:5;
  if caliboptions(jj)==2;
   dkindex(:,jj)=(index+1);
   index=index+1;
  elseif caliboptions(jj)==3;
   ii=(index+1):(index+nbrcam);
   dkindex(:,jj)=ii';
   index=index+nbrcam;
  end;
end;

mut=cell(nbrcam,1);
dt4=[1 0 0]';
dt5=[0 1 0]';
dt6=[0 0 1]';

dK4 =[0 0 1;0 0 0;0 0 0];
dK5 =[0 0 0;0 0 1;0 0 0];
dPdx=sparse(0,0);
for i=1:nbrcam;
  ii=(index+1):(index+6);
  iii=(12*i-11):(12*i);
  dR1=b1*R{i};
  dR2=b2*R{i};
  dR3=b3*R{i};
  KK = [K(i,1) K(i,1)*K(i,3) K(i,4); ...
        0      K(i,1)*K(i,2) K(i,5); ...
        0      0             1 ];  
  
%not:  if i>1, %first camera is fix
  if i>1 | lockcsystem==0,
    dP1 = KK*[dR1 (-dR1*t{i})];
    dP2 = KK*[dR2 (-dR2*t{i})];
    dP3 = KK*[dR3 (-dR3*t{i})];
    dP4 = KK*[zeros(3,3) (-R{i}*dt4)];
    dP5 = KK*[zeros(3,3) (-R{i}*dt5)];
    dP6 = KK*[zeros(3,3) (-R{i}*dt6)];
    dPdx(iii,ii(1))=dP1(:);
    dPdx(iii,ii(2))=dP2(:);
    dPdx(iii,ii(3))=dP3(:);
    if i==2 & lockcsystem>0,
      tmp=[dP4(:),dP5(:),dP6(:)];
      tmp(:,lockcsystem)=[];
      dPdx(iii,ii(4:5))=tmp;
      index=index+5;
    else
      dPdx(iii,ii(4))=dP4(:);
      dPdx(iii,ii(5))=dP5(:);
      dPdx(iii,ii(6))=dP6(:);
      index=index+6;
    end
  end
  
  if dkindex(i,1)>0,
    dK1 =[1 K(i,3) 0;0 K(i,2) 0; 0 0 0];
    dPk1 = dK1*[R{i} (-R{i}*t{i})];
    dPdx(iii,dkindex(i,1)) = dPk1(:);
    if ~isempty(focal),
     if (caliboptions(1)==2 & i==1) | caliboptions(1)==3,
      focalcnt=focalcnt+1;
      focalout(focalcnt,dkindex(i,1))=1/focal(2);
     end
    end
  end;
  if dkindex(i,2)>0,
    dK2 =[0 0 0;0 K(i,1) 0; 0 0 0];
    dPk2 = dK2*[R{i} (-R{i}*t{i})];
    dPdx(iii,dkindex(i,2)) = dPk2(:);
    if ~isempty(ar),
     if (caliboptions(2)==2 & i==1) | caliboptions(2)==3,
      arcnt=arcnt+1;
      arout(arcnt,dkindex(i,2))=1/ar(2);
     end
    end
  end;
  if dkindex(i,3)>0,
    dK3 =[0 K(i,1) 0; 0 0 0; 0 0 0];
    dPk3 = dK3*[R{i} (-R{i}*t{i})];
    dPdx(iii,dkindex(i,3)) = dPk3(:);
  end;
  if dkindex(i,4)>0,
    dPk4 = dK4*[R{i} (-R{i}*t{i})];
    dPdx(iii,dkindex(i,4)) = dPk4(:);
    if ~isempty(pp),
     if (caliboptions(4)==2 & i==1) | caliboptions(4)==3,
      ppcnt=ppcnt+1;
      ppout(ppcnt,dkindex(i,4))=1/pp(3);
     end
    end
  end;
  if dkindex(i,5)>0,
    dPk5 = dK5*[R{i} (-R{i}*t{i})];
    dPdx(iii,dkindex(i,5)) = dPk5(:);
    if ~isempty(pp),
     if (caliboptions(5)==2 & i==1) | caliboptions(5)==3,
      ppcnt=ppcnt+1;
      ppout(ppcnt,dkindex(i,5))=1/pp(3);
     end
    end
  end;
end;





