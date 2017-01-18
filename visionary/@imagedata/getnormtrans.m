function t=getnormtrans(i,option);
% IMAGEDATA/GETNORMTRANS t=getnormtrans(i,option) get normalising image transformation
% INPUT:
%   i - imagedata object
%   option:
%     'affine' - the transformation is restricted to be affine
% OUTPUT:
%   t - 3x3 transformation matrix
% t is chosen such that all coordinates are with [-1,1] in the image.
% If not the image size is available, t is chosen such that the three
% coordinate vectors [x1,x2,...,xn],[y1,y2,...,yn] and [1,1,...,1]
% are orthogonal.
% The transformation should be used to improve numerical conditioning

if nargin<=1,
  option = '';
end

if ~isempty(i.im),

 imsize=size(i.im);
 t =[2/imsize(2)        0  -1
         0    2/imsize(1)  -1
         0         0        1];

else
  %try to use points and conics
  if size(i.conics)>0;
    tmp=[i.points,i.conics(4:6,:)];
  else
    tmp=i.points;
  end
  if size(tmp,2)>2,
    index=find(isfinite(tmp(1,:)));
    tmp=tmp(:,index);
  end
  if size(tmp,2)>2,
    % use points and conics to normalise
    if isempty(strmatch('affine',option))
      [u,s,v]=svd(pflat(tmp));
      t=inv(u*s(:,1:3));
    else
      tmp=pflat(tmp);
      c0=mean(tmp(1:2,:)')';
      mx=max(abs(tmp(1:2,:)-c0*ones(1,size(tmp,2)))');
      d0=diag(1./mx)*sqrt(2);
      t=[[d0,-d0*c0];0,0,1];
    end
  else
    tmp=gethomogeneouslines(i);
    if size(tmp,2)>2,
      index=find(isfinite(tmp(1,:)));
      tmp=tmp(:,index);
    end
    if size(tmp,2)>2,
      % use lines to normalise
      [u,s,v]=svd(pflat(tmp));
      if isempty(strmatch('affine',option))
        t=u*s(:,1:3);
      else
        t=s(:,1:3);
      end
    else
        
        %try again with points and do translation
        if size(i.points,2)==1,
            t=eye(3);
            t(1:2,3)=-i.points(1:2);
        else
          %nothing to normalise with...
          t = eye(3);
        end
        
        
    end
  end
end





