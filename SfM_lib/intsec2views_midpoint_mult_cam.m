function U = intsec2views_midpoint_mult_cam(P,u,cams1,cams2)
% 2 view midpoint-triangulation using midpoint formula.
% P - cameras
% u - image point structure
% cams1 - for each point number of first camera to use
% cams2 - for each point number of second camera to use

%beräkna alla camera centra, vy-riktningar
c1 = NaN(3,u.pointnr);
c2 = NaN(3,u.pointnr);
u1 = NaN(3,u.pointnr);
u2 = NaN(3,u.pointnr);
v1 = NaN(3,u.pointnr);
v2 = NaN(3,u.pointnr);
pp = NaN(3,u.pointnr);
for i = 1:length(P)
    c1(:,cams1 == i) = repmat(-inv(P{i}(:,1:3))*P{i}(:,4),[1 sum(cams1 == i)]);
    c2(:,cams2 == i) = repmat(-inv(P{i}(:,1:3))*P{i}(:,4),[1 sum(cams2 == i)]);
    pp(:,u.index{i}) = u.points{i};
    if sum(~isfinite(pp(1,cams1==i)))
        disp('NaN in cam1');
    end
    if sum(~isfinite(pp(1,cams2==i)))
        disp('NaN in cam2');
    end    
    u1(:,cams1 == i) = pp(:,cams1 == i);
    u2(:,cams2 == i) = pp(:,cams2 == i);
    v1(:,cams1 == i) = inv(P{i}(:,1:3))*u1(:,cams1 == i);
    v2(:,cams2 == i) = inv(P{i}(:,1:3))*u2(:,cams2 == i);
    pp(:,u.index{i}) = NaN;
end
%normalisera vy-riktningar
v1 = v1./repmat(sqrt(sum(v1.^2)),[3,1]);
v2 = v2./repmat(sqrt(sum(v2.^2)),[3,1]);

%normalen för linjerna v1 kryss v2
n = [v1(2,:).*v2(3,:)-v1(3,:).*v2(2,:); ...
     v1(3,:).*v2(1,:)-v1(1,:).*v2(3,:); ...
     v1(1,:).*v2(2,:)-v1(2,:).*v2(1,:)];
n = n./repmat(sqrt(sum(n.^2)),[3 1]);

%kortaste vektor mellan linjerna
d = repmat(sum((c1 - c2).*n),[3 1]).*n;

%Punkt på l2 närmast l1
%planet som innehåller n och l1

npi = [v1(2,:).*n(3,:)-v1(3,:).*n(2,:); ...
     v1(3,:).*n(1,:)-v1(1,:).*n(3,:); ...
     v1(1,:).*n(2,:)-v1(2,:).*n(1,:)];
 
dpi = -sum(c1.*npi);

%Skärning med l2
lambda = (-dpi-sum(c2.*npi))./(sum(npi.*v2));
p2 = c2+repmat(lambda,[3 1]).*v2;

%mittpunkter
U = pextend(p2+d/2);


function y = pextend(x)
y = [x; ones(1,size(x,2))];

