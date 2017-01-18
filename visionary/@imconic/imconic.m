function s=imconic(u,L,n);
% IMCONIC/IMCONIC constructor
% s=imconic(u,L,n);
%

if nargin == 0,
  s.u = [];
  s.L = [];
  s.n = [];
  s = class(s,'imconic');
elseif isa(u,'imconic');
  s=u;
else
  s.u = u;
  s.L = L;
  s.n = n;
  s = class(s,'imconic');
end;
