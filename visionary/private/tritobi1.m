function F=tritobi1(T,im1,im2)

% function F=tritobi1(T,im1,im2)
% Determines the bilinear form between images im1 and im2 from
% the trifocal tensor between images 1, 2 and 3
%
% Input:  T=trifocal tensor between images 1, 2 and 3
%         im1,im2=images for the bilinear form
% Output: F=fundamental matrix between images im1 and im2

ind=sort([im1,im2]);
if ind(1)==2
  im3=1;
elseif ind(2)==3
  im3=2;
else
  im3=3;
end;

for i=1:3
  for j=1:3
    for k=1:3
      ind=[i,j,k];
      ii=ind(im1);
      jj=ind(im2);
      kk=ind(im3);
      TT(ii,jj,kk)=T(i,j,k);
    end;
  end;
end;

F=tritobi(TT);

