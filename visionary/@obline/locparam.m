function obut = locparam(obin,dx);
% OBPOINT/LOCPARAM
%

obut = obin;
dU=obin.dU;
for i=1:4;
 obut.U = obut.U + squeeze(dU(:,:,i))*dx(i);
end;

% Hör borde man normalisera om dvs hitta U som ligger så nära
% U som möjligt i frobeniusnorm och som har ortonormerade kolumner.

