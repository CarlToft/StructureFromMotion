function i=loadimage(i)
%IMAGEDATA/LOADIMAGE loads and stores the image permanently

if ~isempty(i.filename)
 if length(i.filename)>4 & strcmp(i.filename(end-3:end),'.pgm'),
   i.im=readpgm(i.filename);
 else
   i.im = imread(i.filename);
 end
end
