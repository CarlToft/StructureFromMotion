function res=calcresplc(im,s,m,imi,si,mi,offsetres,offsetsm,pp,focal,ar);

antalbilder=size(m,1);
antalfeatures=size(s,2);

res=zeros(offsetres,1);
for ii=1:antalbilder;
 for jj=1:antalfeatures;
  if size(im{ii,jj},1)>0,
   res(imi{ii,jj})=calcres(im{ii,jj},m{ii},s{jj});
  end;
 end;
end;
%focal length residual
if ~isempty(focal),
  for ii=1:antalbilder,
    K=rq(pdp(m{ii}));K=K/K(3,3);
    if (focal(3)==2 & ii==1) | focal(3)==3,
      res(end+1)=(K(1,1)-focal(1))/focal(2);
    end
  end
end
%aspect ratio residual
if ~isempty(ar),
  for ii=1:antalbilder,
    K=rq(pdp(m{ii}));K=K/K(3,3);
    if (ar(3)==2 & ii==1) | ar(3)==3,
      res(end+1)=(K(2,2)/K(1,1)-ar(1))/ar(2);
    end
  end
end
%principal point residual
if ~isempty(pp),
  for ii=1:antalbilder,
    K=rq(pdp(m{ii}));K=K/K(3,3);
    if (pp(4)==2 & ii==1) | pp(4)==3,
      res(end+1)=(K(1,3)-pp(1))/pp(3);
    end
    if (pp(5)==2 & ii==1) | pp(5)==3,
      res(end+1)=(K(2,3)-pp(2))/pp(3);
    end
  end
end
