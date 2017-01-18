function n=addcameras(m,p)
%MOTION/ADDCAMERAS mout=addcameras(m,p) adds p cameras to motion m

n = m;

if iscell(p)
  n.cam = {m.cam{:},p{:}};
else
  n.cam = {m.cam{:},p};
end





