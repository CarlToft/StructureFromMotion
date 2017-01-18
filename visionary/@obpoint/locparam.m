function obut = locparam(obin,dx);
% OBPOINT/LOCPARAM
%

obut = obin;
obut.U = obut.U + obut.dU*dx;
