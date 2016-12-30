function m=changecsystem(m,T,v);
% MOTION/changecsystem m=changecsystem(m,T,v) changes coordinate system
% If T is a 4x4 projective transformation, then P'=P*inv(T) for view v
% If T is a 3x3 projective transformation, then P'=T*P for view v
% If v is not specified all cameras are changed

if nargin<3,
  v=1:length(m.cam);
end

if size(T,1)==4,
  invT=inv(T);
  for i=v,
    m.cam{i} = m.cam{i} * invT;
  end
else
  for i=v,
    m.cam{i} = T*m.cam{i};
  end
end



