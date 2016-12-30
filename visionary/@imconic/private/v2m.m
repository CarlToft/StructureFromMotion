function m=v2m(v);
% V2M m=v2m(v) returns symmetric matrix of conic/quadric vector

if size(v,1)==10,
  m = [ v(1),v(2),v(4),v(7);...
	v(2),v(3),v(5),v(8);...
	v(4),v(5),v(6),v(9);...
	v(7),v(8),v(9),v(10)];

elseif size(v,1)==6,
  m = [ v(1),v(2),v(4);...
	v(2),v(3),v(5);...
	v(4),v(5),v(6)];

else
  error('Wrong vector dimension');
end

