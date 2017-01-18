function i=manpolygon(im,corners,fh);
% i=manline(im,fh) extracts a polygonline in image im.
% The user starts with giving the endpoints of line
% INPUT:
%   im      - image matrix
%   corners - nbr of corners in polygon
%   fh      - (optional) figure handle, if image already plotted
% OUTPUT:
%   i - imagedata object with lines and corner points of polygon

if nargin<=2;
  plot(imagedata(im));
  fh=figure(gcf);
end


figure(fh); hold on;

disp(['Mark the ' num2str(corners) ' corners of the polygon']);
[x,y]=ginput2(corners);
x=[x;x(1)];
y=[y;y(1)];
i=imagedata;
a=1;
troeskel=10;
extracted=zeros(1,corners);

for q=1:corners;
 p1 = [x(q);y(q);1];
 p2 = [x(q+1);y(q+1);1];
 dist = norm(p1-p2);
 antalpunkter=round(dist/5);
 plist = p1*ones(1,antalpunkter) + (p2-p1)*((1:antalpunkter)+2)/(antalpunkter+5);
 n=[0 -1 0;1 0 0;0 0 1]*(p2-p1);
 n=n/norm(n);
 stddevs=zeros(1,antalpunkter);
 normals=zeros(3,antalpunkter);
 edgepoints=zeros(3,antalpunkter);
 k=0;
 for j=1:antalpunkter;
  punkt=plist(:,j);
  sintervall=-6:6;
  [s,sigmas,toppar,mat,matp]=soklangslinje3(im,punkt,n,a,sintervall,troeskel);
  if ~isnan(s),
   k=k+1;
   stddevs(k)=sigmas;
   normals(:,k)=n;
   edgepoints(:,k)=punkt+s*n;
  end;
 end;
 if k<3,
  disp(['Insufficient edge points for a line in polygon']);
 else
  iline=fitline(edgepoints(:,1:k),normals(:,1:k));
  i = i + iline;
  extracted(1,q)=1;
 end;
end;

% calculate corner points
for q=1:corners;
  qprev=mod((q-2),corners)+1;
  if extracted(1,q) & extracted(1,qprev);
    [pn,L]=intersect(i,qprev,q);
    i = addpoints(i,pn,L);
  end
end

%END OF MAIN-FUNCTION


function [pn,L]=intersect(i,line1,line2);
% Intersect two images lines

[l1,L1]=gethomogeneouslines(i,line1);
[l2,L2]=gethomogeneouslines(i,line2);

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

