function i=manconic(im,fh);
% i=manconic(im,fh) extracts a conic in image im.
% The user starts with giving the coordinate of conic centre followed by
%   a point just outside the conic radius.
% INPUT:
%   im - image matrix
%   fh - (optional) figure handle, if image already plotted

if nargin<=1;
  plot(imagedata(im));
  fh=figure(gcf);
end

figure(fh); hold on;
% ask for centre and approx. radius
disp(['Mark center of conic']);
x=ginput2(1);
disp(['Mark approximate radius']);
y=ginput2(1);

dist = norm(x-y);

%extraction
a=1;
troeskel=10;
punkt=[x(1);x(2);1];
antalpunkter=16;
stddevs=zeros(1,antalpunkter);
normals=zeros(3,antalpunkter);
edgepoints=zeros(3,antalpunkter);
k=0;
for theta=(1:antalpunkter)*2*pi/antalpunkter;
 n=[cos(theta);sin(theta);0];
 sintervall=0:round(dist);
 [s,sigmas,toppar]=soklangslinje3(im,punkt,n,a,sintervall,troeskel);
 if ~isnan(s),
   k=k+1;
   stddevs(k)=sigmas;
   normals(:,k)=n;
   edgepoints(:,k)=punkt+s*n;
 end
end;

if k<6,
  disp(['Could not extract sufficient number of edgepoints']);
  i=imagedata;
else
  %fit conic with normal
  i=fitconic(edgepoints(:,1:k),normals(:,1:k));
end

%END OF MAIN FUNCTION
