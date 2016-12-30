function p=getcameras(m,index)
%MOTION/GETCAMERAS getcameras(m,index) returns m's cameras with index

if nargin==1,
  p=m.cam;
else
  if length(index)==1,
    p=m.cam{index};
  else
    p={m.cam{index}};
  end
end



