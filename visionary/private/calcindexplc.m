function [mi,si,imi,offsetmotion,offsetsm,offsetres]=calcindexplc(im,s,m,mode,zplane);
%
% calculates index structures needed in bundle adjustment.
% INPUT
%   im - cell matrix of @im-feature objects
%   s  - cell vector of @ob-feature objects
%   m  - ell vector of @pmatrix objects
% OUTPUT
% se själv


antalbilder=size(m,1);
antalfeatures=size(s,2);

mi=cell(antalbilder,1);
si=cell(1,antalfeatures);
imi=cell(antalbilder,antalfeatures);

offset=0;
if mode~=4 & mode~=7, %only structure
  for ii=1:antalbilder;
    sx = sizedx(m{ii});
    mi{ii} = (offset+1):(offset+sx);
    offset=offset+sx;
  end;
elseif mode==7, %known rotation
  for ii=1:antalbilder;
%    sx = sizedx(m{ii});
    sx = 3;
    mi{ii} = (offset+1):(offset+sx);
    offset=offset+sx;
  end;
end    

offsetmotion=offset;

if mode~=5 & mode ~=6, %only motion
  for ii=1:antalfeatures;
    sx = sizedx(s{ii});
    if ~isempty(find(ii==zplane))>0,
        sx=sx-1;
    end
    si{ii} = (offset+1):(offset+sx);
    offset=offset+sx;
  end;
end

offsetsm=offset;

offset=0;
for ii=1:antalbilder;
 for jj=1:antalfeatures;
  if size(im{ii,jj},1)>0,
   sr = sizeres(im{ii,jj});
   imi{ii,jj} = (offset+1):(offset+sr);
   offset=offset+sr;
  end;
 end;
end;
offsetres=offset;

