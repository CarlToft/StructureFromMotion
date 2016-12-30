function [r,drdP,drdU]=dres(im,proj,ob);
% IMPOINT/DRES function [r,drdP,drdU]=dres(im,proj,ob);
% calculates weighted residuals and their derivatives with
% respect to changes in structure and motion

[U,dU]=udu(ob);

[P,dP]=pdp(proj);

up = P*U;
dupdP = squeeze(dP(:,1,:)*U(1) + dP(:,2,:)*U(2) + dP(:,3,:)*U(3) + dP(:,4,:)*U(4));
dupdU = P*dU;

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
