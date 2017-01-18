function s=pmatrix(P);

if nargin == 0,
  s.P = [];
  s.dP = [];
  s = class(s,'pmatrix');
elseif isa(P,'pmatrix');
  s = P;
else
  s.P = P;
  s.dP = reshape(eye(12,12),3,4,12);
  s = class(s,'pmatrix');
end;
