function apoint=intersect(line1,line2);

l1=line1.u;
l2=line2.u;
L1=line1.L;
L2=line2.L;
C1=pinv(L1'*L1);
C2=pinv(L2'*L2);

T1=[0 -l1(3) l1(2);l1(3) 0 -l1(1);-l1(2) l1(1) 0];
T2=[0 -l2(3) l2(2);l2(3) 0 -l2(1);-l2(2) l2(1) 0];

p1 = cross(l1,l2);
Cp = T1*C2*T1' + T2*C1*T2';

n=[0;0;1];
pn= p1/(p1'*n);
dpndp1 = (eye(3)/(p1'*n)) - (p1*n')/(p1'*n)^2;
Cpn = dpndp1*Cp*dpndp1';
[U,S,V]=svd(Cpn);
L=sqrtm(inv(S(1:2,1:2)))*U(:,1:2)';
apoint = impoint(pn,L,n);
