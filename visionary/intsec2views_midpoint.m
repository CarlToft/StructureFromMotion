function U = intsec2views_midpoint(P1,P2,u1,u2)
%beräkna linjerna
v1 = inv(P1(:,1:3))*u1;
v1 = v1./repmat(sqrt(sum(v1.^2)),[3,1]);
c1 = -inv(P1(:,1:3))*P1(:,4);

v2 = inv(P2(:,1:3))*u2;
v2 = v2./repmat(sqrt(sum(v2.^2)),[3,1]);
c2 = -inv(P2(:,1:3))*P2(:,4);

%normalen för linjerna v1 kryss v2
n = [v1(2,:).*v2(3,:)-v1(3,:).*v2(2,:); ...
     v1(3,:).*v2(1,:)-v1(1,:).*v2(3,:); ...
     v1(1,:).*v2(2,:)-v1(2,:).*v2(1,:)];
n = n./repmat(sqrt(sum(n.^2)),[3 1]);

%kortaste vektor mellan linjerna
d = repmat((c1 - c2)'*n,[3 1]).*n;

%Punkt på l2 närmast l1
%planet som innehåller n och l1

npi = [v1(2,:).*n(3,:)-v1(3,:).*n(2,:); ...
     v1(3,:).*n(1,:)-v1(1,:).*n(3,:); ...
     v1(1,:).*n(2,:)-v1(2,:).*n(1,:)];
 
dpi = -c1'*npi;

%Skärning med l2
lambda = (-dpi-c2'*npi)./(sum(npi.*v2));
p2 = repmat(c2,[1 size(v2,2)])+repmat(lambda,[3 1]).*v2;

%mittpunkter
U = [p2+d/2; ones(1,size(p2,2))];


