function [r]=calcres(im,proj,ob);
% IMLINE/CALCRES function [r]=res(im,proj,ob);
% calculates weighted residuals

[U,dU]=udu(ob);

[P,dP]=pdp(proj);

upl = pflat(P*U);


r = upl'*im.u/im.stddevs;

%up = pflat(cross(upl(:,1),upl(:,2)));
%L=im.L;
%n=im.n;
%un= up/(up'*n);
%r = L*(um-un);
