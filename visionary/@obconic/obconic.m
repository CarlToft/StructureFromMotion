function s=obconic(U);
% OBCONIC/OBCONIC constructor
%

if nargin == 0,
  s.U = [];
  s.dU = [];
  s = class(s,'obconic');
elseif isa(U,'obconic');
  s = U;
else
  U = U/norm(U,'fro');
  [uu,ss,vv]=svd(m2v(U));
  dU=zeros(4,4,9);
  for i=1:9; dU(:,:,i)=v2m(uu(:,i+1)); end;
  s.U = U;
  s.dU = dU;
  s = class(s,'obconic');
end;
