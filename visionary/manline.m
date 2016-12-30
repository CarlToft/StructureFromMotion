function i=manline(im,fh);
% i=manline(im,fh) extracts a line in image im.
% The user starts with giving the endpoints of line
% INPUT:
%   im - image matrix
%   fh - (optional) figure handle, if image already plotted

if nargin<=1;
  plot(imagedata(im));
  fh=figure(gcf);
end

figure(fh); hold on;
% ask for centre and approx. radius
disp(['Mark endpoints of line']);
[x,y]=ginput2(2);


p1 = [x(1);y(1);1];
p2 = [x(2);y(2);1];
dist = norm(p1-p2);
antalpunkter=round(dist/5);
plist = p1*ones(1,antalpunkter)+(p2-p1)*((1:antalpunkter)+2)/(antalpunkter+5);
n=[0 -1 0;1 0 0;0 0 1]*(p2-p1);
n=n/norm(n);
stddevs=zeros(1,antalpunkter);
normals=zeros(3,antalpunkter);
edgepoints=zeros(3,antalpunkter);
a=1;
troeskel=10;
k=0;
for j=1:antalpunkter;
  punkt=plist(:,j);
  sintervall=-6:6;
  [s,sigmas,toppar]=soklangslinje3(im,punkt,n,a,sintervall,troeskel);
  if ~isnan(s),
   k=k+1;
   stddevs(k)=sigmas;
   normals(:,k)=n;
   edgepoints(:,k)=punkt+s*n;
  end
end;

if k<3,
  disp(['Could not extract sufficient number of edgepoints']);
  i=imagedata;
else
  %fit line with normal
  i=fitline(edgepoints(:,1:k),normals(:,1:k));
end

%END OF MAIN FUNCTION
