function intensities = calcintensities(basep,bild,H);

[intensities, intengrad]=measure(bild,pflat(H*basep.points),basep.a);
