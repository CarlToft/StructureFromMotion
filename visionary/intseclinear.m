function s=intseclinear(m,imseq,option,soriginal);
% INTSECLINEAR s=intseclinear(m,imseq,option,soriginal) calculates 3D features
% with linear methods
% INPUT:
%   m - motion object
%   imseq - cell array of imagedata objects
%   soriginal - If specified, only "missing" features in this object are intersected
%   option:
%     'nocoordtrans' - No coordinate transformation
%     'conic' - Quadrics are intersected as conics
% OUTPUT:
%   s - structure object containing reconstructed 3D features
% The results may be inaccurate.


if nargin<=2
  option = [];
end

if nargin<=3;
  soriginal=[];
end

s = intsecpoints(m,imseq,soriginal,option);
s = s + intseclines(m,imseq,soriginal);
s = s + intsecconics(m,imseq,soriginal,option);


