function [r]=calcres(im,proj,ob,imnbr);
% IMPOINT/CALCRES function [r]=res(im,proj,ob);
% calculates weighted residuals

%[U,dU]=udu(ob);
U=udu(ob);
%[P,dP]=pdp(proj);
P=pdp(proj);

up = P*U;

n=im.n;
um=im.u;
L=im.L;

un= up/(up'*n);

r = L*(um-un);
