function A=skew(a);
% A=skew(a) - returns skew matrix y such that a x v = Av for any v
%

A=[0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0];
