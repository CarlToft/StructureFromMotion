function [r]=plotres(im,proj,ob);
% IMCONIC/CALCRES function [r]=res(im,proj,ob);
% calculates weighted residuals

[U,dU]=udu(ob);

[P,dP]=pdp(proj);

up = m2v(P*U*P');

n=im.n;
um=im.u;
L=im.L;

un= up/(up'*n);
ritac(un,'k--');

r = L*(um-un);
