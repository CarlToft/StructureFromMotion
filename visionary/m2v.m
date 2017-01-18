function v=m2v(m);
% M2V v=m2v(m) returns vector of symmetric conic/quadric matrix

if size(m,1)==4,
   v=[m(1,1),m(1,2),m(2,2),m(1,3),m(2,3),m(3,3),...
	m(1,4),m(2,4),m(3,4),m(4,4)]';
elseif size(m,1)==3,
   v=[m(1,1),m(1,2),m(2,2),m(1,3),m(2,3),m(3,3)]';
else
  error('Wrong matrix-dimension');
end

      
