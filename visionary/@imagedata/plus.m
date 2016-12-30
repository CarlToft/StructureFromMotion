function i=plus(i1,i2)
%IMAGEDATA/PLUS i=plus(i1,i2) adds i2's features to i1

i = imagedata(i1);
i2 = imagedata(i2);

[p,pc]=getpoints(i2);
i.points = [i.points,p];
i.pointcov = [i.pointcov,pc];

[p,pc]=getlines(i2);
if isempty(i.lines),
  i.lines = [p];
  i.linecov = [pc];
else
  i.lines = [i.lines,p];
  i.linecov = [i.linecov,pc];
end
  
[p,pc]=getconics(i2);
i.conics = [i.conics,p];
i.coniccov = [i.coniccov,pc];

