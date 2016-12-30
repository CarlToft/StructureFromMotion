function obut = setparam(obin,xparam);
% OBPOINT/SETPARAM
%

obut = obin;
currparam = getparam(obin);

obut.U = obut.U + obut.dU*(xparam-currparam);
