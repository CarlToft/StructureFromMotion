function i = imagedata(im,points,lines,conics)
%IMAGEDATA class constructor
%  i=IMAGEDATA(im,points,lines,conics) creates an IMAGEDATA object where
%   im: Filename of image or image matrix
%   points : 3xn (or 2xn) matrix where columns are homogeneous
%            (or cartesian) point coordinates
%   lines  : 3xm matrix where columns are homogeneous line coordinates
%            (or 6xm matrix where columns are the two homogeneous end points
%             of the line)
%   conics : 6xp matrix where columns are dual conic coordinates
%
%  All of the input arguments are optional
%
% Methods:
%
%  addconics
%  addlines
%  addpoints
%  changecsystem
%  clearconics
%  clearimage
%  clearlines
%  clearpoints
%  display
%  getconics
%  getfilename
%  getimage
%  getlines
%  getnormtrans
%  getpoints
%  imagedata
%  loadimage
%  plot
%  plus
%  size


if nargin == 1 & isa(im,'imagedata'),
 i=im;
else
 i.filename = [];
 i.im = [];
 i.pointcov={};
 i.linecov={};
 i.coniccov={};

 if nargin >= 1 & ~isempty(im),
   if ischar(im),
     i.filename = im;
   else
     i.im = im;
   end
 end
 if nargin >= 2 & ~isempty(points),
   if size(points,1)==2.
    i.points=[points;ones(1,size(points,2))];
   else
    i.points = points;
   end
 else
   i.points = zeros(3,0);
 end
 if nargin >= 3 & ~isempty(lines),
   i.lines = lines;
 else
   i.lines = zeros(3,0);
 end
 if nargin >= 4 & ~isempty(conics),
   i.conics = conics;
 else
   i.conics = zeros(6,0);
 end

 i = class(i,'imagedata');

end


 
