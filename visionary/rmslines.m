function tmp=rmslines(s,m,imseq,v);
% RMSLINES rmslines(s,m,imseq,v) reproject structure s
%  with motion m and calculate Root Mean Square (RMS) error for lines in images
%  Residuals are distance of project 3D end points to line.
% INPUT:
%  s - structure
%  m - motion
%  imseq - cell list of imagedata objects
%  v - view . If not specified, the whole sequence is displayed
%
% See also REPROJECT

if nargin <=3,
  vs = 1:size(m);
else
  vs = v;
end
nbrvs=length(vs);
rms=zeros(1,nbrvs);
nbr=zeros(1,nbrvs);
for i=vs;
  imreproj=project(s,m,i);
  lll=gethomogeneouslines(imseq{i});
  reendpts=getlines(imreproj);
  reendpts(1:3,:)=pflat(reendpts(1:3,:));
  reendpts(4:6,:)=pflat(reendpts(4:6,:));
  lindex=find(~isnan(lll(1,:)) & ~isnan(reendpts(1,:)));
  nbrlll=length(lindex);
  err=zeros(1,2*length(lindex));
  for jj=1:length(lindex);
    err(2*jj-1:2*jj)=lll(:,lindex(jj))'/norm(lll(1:2,lindex(jj)))*[reendpts(1:3,lindex(jj)),reendpts(4:6,lindex(jj))];
  end
  
  rms(i)=sum(err.^2);
  tmp=sqrt(rms(i)/(2*nbrlll));
  if nargout==0,
    disp(['View ',num2str(i),', RMS ',num2str(tmp),' --- Number of lines: ',num2str(nbrlll)]);
  end
  nbr(i)=nbrlll;
end
tmp=sqrt(sum(rms)/(2*sum(nbr)));
if nargout ==0,
  disp(['Total RMS ',num2str(tmp),' --- Number of lines: ',num2str(sum(nbr))]);
end


