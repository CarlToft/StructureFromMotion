function [inten, intengrad]=measure(image,points,a);

inten=0;

for i=1:size(points,2)
  [In, Igrad] = measurepoint(image,points(1:2,i),a);
  inten(1,i) = In;
  intengrad(1:2,i) = Igrad;
end;

