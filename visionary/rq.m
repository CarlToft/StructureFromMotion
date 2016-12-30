function [r,q]=rq(a)
% RQ [r,q]=rq(a) factorises a such a=rq where r upper tri. and q unit matrix
% If a is not square, then q is equal q=[q1 q2] where q1 is unit matrix

[m,n]=size(a);
e=eye(m);
p=e(:,[m:-1:1]);
[q0,r0]=qr(p*a(:,1:m)'*p);

r=p*r0'*p;
q=p*q0'*p;

fix=diag(sign(diag(r)));
r=r*fix;
q=fix*q;

if n>m
  q=[q, inv(r)*a(:,m+1:n)];
end
