function s=obline(U);
% OBLINE/OBLINE constructor
%

if nargin == 0,
  s.U = [];
  s.dU = [];
  s = class(s,'obline');
elseif isa(U,'obline');
  s = U;
else
%  U=psphere(U); %not valid with end point norm
%  du=null(U');
%  du1=du(:,1);
%  du2=du(:,2);
  U = pflat(U);
  dir= U(:,2)-U(:,1);
  [du,ds,dv]=svd([dir [0;0;0;1]]);
  du1=du(:,3);
  du2=du(:,4);
  
  s.U = U;
  s.dU = cat(3,[du1 zeros(4,1)],[du2 zeros(4,1)],[zeros(4,1) du1],[zeros(4,1) du2]);
  s = class(s,'obline');
end;
