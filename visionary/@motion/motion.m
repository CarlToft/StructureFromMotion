function m = motion(p)
%MOTION class constructor
%  m=motion(p) creates an MOTION object where
%   p: is a cell array of 3x4 camera matrices or single 3x4 camera matrix

if nargin == 1 & isa(p,'motion'),
 m=p;
else

 if nargin == 0,
  m.cam = {};
 else
  if iscell(p),
    m.cam = p;
  else
      if size(p,1)==3 & size(p,2)==4,
          m.cam = {p};
      else
          m.cam = {};
      end
  end
 end
 m = class(m,'motion');

end
