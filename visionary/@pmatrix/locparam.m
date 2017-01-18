function obut = locparam(obin,dx);
% PMATRIX/LOCPARAM
%
obut = obin;
dP=obut.dP;
for i=1:size(dx,1);
 obut.P = obut.P + squeeze(dP(:,:,i))*dx(i);
end;

