function p=findproj(data1,data2);
% p=FINDPROJ(data1,data2) - finds projective transformation 'p' using points
% "p(data1) = data2"
% data1, data2 : structure or imagedata objects
% See also: changecsystem

b1 = getpoints(data1);
b2 = getpoints(data2);

M1=[];M2=[];
[m,n]=size(b1);
for i=1:m,
 M1=[M1; zeros(n,m*(i-1)) b1' zeros(n,m*(m-i))];
 M2=[M2; diag(b2(i,:))];
end;
M=[M1 (-M2)];
[Q,R]=qr(M);
A=R(1:(n+m*m-1),1:(n+m*m-1));
b=R(1:(n+m*m-1),n+m*m);
x=[inv(A)*b;-1];
p=zeros(m,m);
p(:)=x(1:m*m);
p=p';
