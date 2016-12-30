function obut = locparam(obin,dx);
% OBPOINT/LOCPARAM
%

obut = obin;
dU=obin.dU;
for i=1:9;
 obut.U = obut.U + squeeze(dU(:,:,i))*dx(i);
end;

% H�r borde man normalisera om dvs hitta U som ligger s� n�ra
% U som m�jligt i frobeniusnorm och som har ortonormerade kolumner.

