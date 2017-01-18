function s=obpoint(U);

if nargin == 0,
  s.U = [];
  s.dU = [];
  s = class(s,'obpoint');
elseif isa(U,'obpoint');
  s = U;
else
%  s.U = pflat(U);  %not necessary - causes numerical instabilities
  s.U = U;
  s.dU = [1 0 0;0 1 0;0 0 1;0 0 0];
  s = class(s,'obpoint');
end;
