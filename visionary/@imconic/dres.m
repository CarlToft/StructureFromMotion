function [r,drdP,drdU]=dres(im,proj,ob);
% IMCONIC/DRES function [r,drdP,drdU]=dres(im,proj,ob);
% calculates weighted residuals and their derivatives with
% respect to changes in structure and motion

[U,dU]=udu(ob);

[P,dP]=pdp(proj);

up = m2v(P*U*P');

dupdU=zeros(6,sizedx(ob));
for i=1:sizedx(ob);
  dupdU(:,i)=m2v(P*squeeze(dU(:,:,i))*P');
end;
dupdP=zeros(6,sizedx(proj));
for i=1:sizedx(proj);
  dupdP(:,i)=m2v( squeeze(dP(:,:,i))*U*P' + P*U*squeeze(dP(:,:,i))' );
end;

n=im.n;
um=im.u;
L=im.L;

un= up/(up'*n);
dundup = (eye(size(up,1))/(up'*n) - up*n'/(up'*n)^2);
dundP = dundup*dupdP;
dundU = dundup*dupdU;

r = L*(um-un);
drdP = -L*dundP;
drdU = -L*dundU;

%keyboard;
