function m=plus(m1,m2)
%MOTION/PLUS m=plus(m1,m2) adds m2's cameras to m1

m = motion(m1);
m2 = motion(m2);

m.cam = [m.cam, getcameras(m2)];


