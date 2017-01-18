function [i,crf]=harrispro(im,option);
% i=HARRISPRO(im,option) - creates an IMAGEDATA object with corners
% using HARRIS, plus auto-correlation in neighbourhood
% im: image matrix
% option: (cell of strings)
% - 'threshold=X', sets corner response threshold to X
% - 'nbrcorners=X', gets the X strongest corners in the image
% - 'noimage', image is not included in IMAGEDATA object i
% Output:
%  i - IMAGEDATA object
%  crf - corner response function (vector)
% Default is 'threshold=1'
% See also: harris.m

ignore=3;
autothreshold=0.2;
patchsize=5;
searchsize=50;

if nargin<2,
  option=[];
end

i1=harris(im,option);

n=2*patchsize+1;
imconv=conv2(im.^2,ones(n),'valid')-conv2(im,ones(n),'valid').^2/n^2;    
imconv(find(imconv==0))=1;

i=imagedata;
crf=[];

for jj=1:size(i1,1);

  pos=round(getpoints(i1,jj));
  patch=im(pos(2)+(-patchsize:patchsize),pos(1)+(-patchsize:patchsize));
  patch = flipud(fliplr(patch));
  patch=patch-mean(patch(:));
  patch2=sum(patch(:).^2);if patch2==0,patch2=1;end;

  index1=max(1,pos(2)-searchsize):min(size(im,1),pos(2)+searchsize);
  index2=max(1,pos(1)-searchsize):min(size(im,2),pos(1)+searchsize);
  searchpatch=im(index1,index2);


  searchsquared=imconv(index1(1+patchsize:end-patchsize)-patchsize,...
		 index2(1+patchsize:end-patchsize)-patchsize);


  s_patch=conv2(searchpatch,patch,'valid');

  res=(searchsquared+patch2-2*s_patch)./sqrt(searchsquared)/sqrt(patch2);

  res(pos(2)-index1(1)+1-patchsize+(-ignore:ignore),...
      pos(1)-index2(1)+1-patchsize+(-ignore:ignore))=Inf;

  residual=min(res(:));

  if residual>autothreshold,
    i=addpoints(i,pos);
    crf=[crf;residual];
  end
end

