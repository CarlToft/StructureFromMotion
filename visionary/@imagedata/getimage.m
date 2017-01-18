function im=getimage(i)
%IMAGEDATA/GETIMAGE im=getimage(i) returns i's image matrix.

if ~isempty(i.filename) & isempty(i.im),
  if length(i.filename)>4 & strcmp(i.filename(end-3:end),'.pgm'),
    im=readpgm(i.filename);
  else
    im = imread(i.filename);
  end
else
  im=i.im;
end



