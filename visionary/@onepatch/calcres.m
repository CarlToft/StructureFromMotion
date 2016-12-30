function [res,intm,intp] = calcres(basep,bild,H,dhdx);

[intensities, intengrad]=measure(bild,pflat(H*basep.points),basep.a);

%dHdx = homographyderiv(H,dhdx);
pph = H*basep.points;
pp = pflat(pph);

[intensities, intengrad]=measure(bild,pp,basep.a);

res = (intensities - basep.intensities)';
intm = intensities';
intp = basep.intensities';
