function s=impoint(u,L,n);

if nargin == 0,
  s.u = [];
  s.L = [];
  s.n = [];
  s = class(s,'impoint');
elseif isa(u,'impoint');
  s = u;
else
  s.u = u;
  s.L = L;
  s.n = n;
  s = class(s,'impoint');
end;
