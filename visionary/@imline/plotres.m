function [r]=plotres(im,proj,ob);
% IMLINE/PLOTRES function [r]=res(im,proj,ob);
% calculates weighted residuals

%keyboard;
[U,dU]=udu(ob);

[P,dP]=pdp(proj);

upl = P*U;

up = cross(upl(:,1),upl(:,2));

n=im.n;
um=im.u;
L=im.L;

un= up/(up'*n);
rital(un,'k--');

r = L*(um-un);
