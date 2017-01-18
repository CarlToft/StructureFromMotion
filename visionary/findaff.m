function a=findaff(data1,data2);
% a=FINDAFF(data1,data2) - finds affine transformation 'a' using points
% "a(data1) = data2"
% data1, data2 : structure or imagedata objects
% See also: changecsystem

x = pflat(getpoints(data1));
y = pflat(getpoints(data2));

m = size(x,1)-1;
x = x(1:m,:);
y = y(1:m,:);

xmean=mean(x')';
ymean=mean(y')';

dx = x-xmean*ones(1,size(x,2));
dy = y-ymean*ones(1,size(y,2));

A = dy/dx;

a = [A,ymean-A*xmean;zeros(size(xmean))',1];

