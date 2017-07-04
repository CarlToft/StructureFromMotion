

im=imread('blueberry.tif');
%figure(1);clf;image(im);

path = 'results';
[kp, files] = readKPandDesc(path)


locs = kp.locs;
figure(1);clf;
plot(imagedata(im,locs));
