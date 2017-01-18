function [r]=plotres(im,proj,ob);
% IMPOINT/PLOTRES function [r]=res(im,proj,ob);
% calculates weighted residuals

[U,dU]=udu(ob);

[P,dP]=pdp(proj);

up = P*U;

n=im.n;
um=im.u;
L=im.L;

un= up/(up'*n);
rita(pflat(un),'y*');
rita(pflat(um),'r*');

r = L*(um-un);
