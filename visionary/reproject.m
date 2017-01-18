function reproject(s,m,imseq,v,option);
% REPROJECT reproject(s,m,imseq,v,option) reprojects structure s
%  with motion m and compares to imseq
%  v - view . If not specified, the whole sequence is displayed
%  option:
%    'numbered' - points are numbered
%
% See also RMSPOINTS

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
if nargin<5,
  option='';
end


for i=vs;
  imreproj=project(s,m,i);
  clf;
  if isa(imseq{i},'imagedata'),
      imdata=imseq{i};
  else
      imdata=imagedata([],imseq{i});
  end
  plot(imdata);
  [p,l,c]=getcommonfeatures({imdata,imreproj});
  imreproj=imagedata([],getpoints(imreproj,p),getlines(imreproj,l),getconics(imreproj,c));
  hold on; zoom on;
  plot(imreproj,'ro','r-.','r-.');
  if strcmp(option,'numbered'),
    pts=pflat(getpoints(imreproj));
    ax=axis;
    for pp=1:size(pts,2);dev=0.01*ax(2);
      h=text(pts(1,pp)-dev,pts(2,pp)-dev,num2str(p(pp)));
      set(h,'Fontsize',11);set(h,'Color','b');
    end
  end
  if i ~= vs(length(vs)),
    disp(['View ' num2str(i) '. Press space.']);
    pause;
  else
    disp(['View ' num2str(i) '.']);
  end
  hold off;
end
