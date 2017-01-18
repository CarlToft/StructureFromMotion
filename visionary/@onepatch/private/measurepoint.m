function [inten, intengrad]=measurepoint(bild,punkt,a);

if a<1
  error('Implementation only supports a-values larger than 1');
end;

x0=punkt(1);
y0=punkt(2);
[m,n]=size(bild);
NN=round(a*3);


cutx=max(round(x0)-NN,1):min(round(x0)+NN,n);
cuty=max(round(y0)-NN,1):min(round(y0)+NN,m);
cutbild=bild(cuty,cutx);

[x,y]=meshgrid(cutx,cuty);

filter=exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^2*pi);
filtx=2*(x-x0).*exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^4*pi);
filty=2*(y-y0).*exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^4*pi);
intengrad(1,1)=sum(sum(filtx.*cutbild));
intengrad(2,1)=sum(sum(filty.*cutbild));
inten=sum(sum(filter.*cutbild));

