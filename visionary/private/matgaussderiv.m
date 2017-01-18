function [matning,komplext]=matgaussderiv(bild,punkt,theta,a);

x0=punkt(1);
y0=punkt(2);
[m,n]=size(bild);
NN=round(a*3);
cutx=max(round(x0)-NN,1):min(round(x0)+NN,n);
cuty=max(round(y0)-NN,1):min(round(y0)+NN,m);
cutbild=bild(cuty,cutx);
[x,y]=meshgrid(cutx,cuty);
filtx=2*(x-x0).*exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^4*pi);
filty=2*(y-y0).*exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^4*pi);
matx=sum(sum(filtx.*cutbild));
maty=sum(sum(filty.*cutbild));
matning=cos(theta)*matx+sin(theta)*maty;
komplext=matx+sqrt(-1)*maty;


