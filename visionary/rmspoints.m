function tmp=rmspoints(s,m,imseq,v);
% RMSPOINTS rmspoints(s,m,imseq,v) reproject structure s
%  with motion m and calculate Root Mean Square (RMS) error for points in images
% INPUT:
%  s - structure
%  m - motion
%  imseq - cell list of imagedata objects
%  v - view . If not specified, the whole sequence is displayed
%
% See also REPROJECT

if ~isa(s,'structure'),
    s=structure(s);
end
if ~isa(m,'motion'),
    m=motion(m);
end

if nargin <=3 | isempty(v),
  vs = 1:size(m);
else
  vs = v;
end
nbrvs=length(vs);
rms=zeros(1,nbrvs);
nbr=zeros(1,nbrvs);
res=[];
for i=vs;
  imreproj=project(s,m,i);
  if isa(imseq{i},'imagedata'),
      imdata=imseq{i};
  else
      imdata=imagedata([],imseq{i});
  end

  pts=pflat(imdata);
  repts=pflat(imreproj);
  pindex=find(~isnan(pts(1,:)) & ~isnan(repts(1,:)));
  nbrpts=length(pindex);
  err=repts(1:2,pindex)-pts(1:2,pindex);
  res=[res;err(:)];

  rms(i)=sum(sum(err.^2));
  tmp=sqrt(rms(i)/(2*nbrpts));
  if nargout==0,
    disp(['View ',num2str(i),', RMS ',num2str(tmp),' --- Number of points: ',num2str(nbrpts)]);
  end
  nbr(i)=nbrpts;
end
tmp=sqrt(sum(rms)/(2*sum(nbr)));
if nargout ==0,
  disp(['Total RMS ',num2str(tmp),' --- Number of points: ',num2str(sum(nbr))]);
  clf;plot(res,'.');
end
