function [sut,mut]=dolocparamplc(s,m,si,mi,dx,zplane);

antalbilder=size(m,1);
antalfeatures=size(s,2);

sut=cell(1,antalfeatures);
mut=cell(antalbilder,1);

for ii=1:antalbilder;
  if ~isempty(mi{ii}),
      if length(mi{ii})==3, %knownrotation
          dxtmp=zeros(12,1);
          dxtmp(10:12)=dx(mi{ii});
          mut{ii}=locparam(m{ii},dxtmp);
      else
          mut{ii}=locparam(m{ii},dx(mi{ii}));
      end
  else
    mut{ii}=m{ii};
  end
end;

for ii=1:antalfeatures;
  if ~isempty(si{ii}),
    if ~isempty(find(ii==zplane))>0,
      sut{ii}=locparam(s{ii},[dx(si{ii});0]);
    else
      sut{ii}=locparam(s{ii},dx(si{ii}));
    end
        
      
  else
    sut{ii}=s{ii};
  end
end;

