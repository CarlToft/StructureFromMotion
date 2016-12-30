function  [mut] = aclocparam(K,R,t,dy,caliboptions,lockcsystem);

 b1=[0 1 0;-1 0 0;0 0 0];
 b2=[0 0 -1;0 0 0;1 0 0];
 b3=[0 0 0;0 0 1;0 -1 0];

nbrcam=size(K,1);
index=0;
for jj=1:5;
  if caliboptions(jj)==2;
   ii=index+1;
   K(:,jj)=K(:,jj)+dy(ii);
   index=index+1;
  elseif caliboptions(jj)==3;
   ii=(index+1):(index+nbrcam);
   K(:,jj)=K(:,jj)+dy(ii);
   index=index+nbrcam;
  end;
end;

mut=cell(nbrcam,1);
for i=1:nbrcam;
  if i>1 | lockcsystem==0, %first camera is fix
    if i==2 & lockcsystem>0,
      ii=(index+1):(index+5);
      dxx=dy(ii);
      indt=1:3;indt(lockcsystem)=[];
      dt=zeros(3,1);
      dt(indt)=dxx(4:5);
      t{i}=t{i} + dt;
      index=index+5;
    else
      ii=(index+1):(index+6);
      dxx=dy(ii);
      t{i}=t{i} + dxx(4:6);
      index=index+6;
    end
    R{i}=expm(dxx(1)*b1+dxx(2)*b2+dxx(3)*b3)*R{i};
  end
  KK = [K(i,1) K(i,1)*K(i,3) K(i,4); ...
        0      K(i,1)*K(i,2) K(i,5); ...
        0      0             1 ];
  mut{i}=pmatrix(KK*[R{i} (-R{i}*t{i})]);
end;
