function [res,resgrad]=calcresgradplc(im,s,m,imi,si,mi,offsetres,offsetsm,pp,focal,ar,zplane);

antalbilder=size(m,1);
antalfeatures=size(s,2);

res=zeros(offsetres,1);
resgrad=sparse(offsetres,offsetsm);
for ii=1:antalbilder;
 for jj=1:antalfeatures;
%  [r,drdP,drdU]=dres(im{ii,jj},m{ii},s{jj});
  if size(im{ii,jj},1)>0,
   [res(imi{ii,jj}),motgrad,strgrad]=dres(im{ii,jj},m{ii},s{jj});
   if ~isempty(mi{ii}),
       if length(mi{ii})==3, %knownrotation case
           resgrad(imi{ii,jj},mi{ii})=motgrad(:,10:12);
       else
           resgrad(imi{ii,jj},mi{ii})=motgrad;
       end
   end
   if ~isempty(si{jj}),
     if ~isempty(find(jj==zplane))>0,
         strgrad(:,3)=[];
     end
     resgrad(imi{ii,jj},si{jj})=strgrad;
   end
  end;
 end;
end;

%focal length residual (the gradpart is in accalcgrad.m)
if ~isempty(focal),
  for ii=1:antalbilder,
    K=rq(pdp(m{ii}));K=K/K(3,3);
    if (focal(3)==2 & ii==1) | focal(3)==3,
      res(end+1)=(K(1,1)-focal(1))/focal(2);
    end
  end
end
%aspect ratio residual (the gradpart is in accalcgrad.m)
if ~isempty(ar),
  for ii=1:antalbilder,
    K=rq(pdp(m{ii}));K=K/K(3,3);
    if (ar(3)==2 & ii==1) | ar(3)==3,
      res(end+1)=(K(2,2)/K(1,1)-ar(1))/ar(2);
    end
  end
end
%principal point residual (the gradpart is in accalcgrad.m)
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
