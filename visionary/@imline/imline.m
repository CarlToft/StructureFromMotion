function s=imline(u,stddevs);
% IMLINE/IMLINE constructor
% s=imline(u,stddevs);
%


%old: function s=imline(u,L,n);

if nargin == 0,
  s.u = [];
  s.stddevs = [];
%  s.L = [];
%  s.n = [];
  s = class(s,'imline');
elseif isa(u,'imline');
  s = u;
else
  s.u = u;
  s.stddevs = stddevs;
  s = class(s,'imline');
%  s.L = L;
%  s.n = n;
end;
