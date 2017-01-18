function plot(i,p,l,c)
%IMAGEDATA/PLOT  Plot imagedata.
%  PLOT(i,p,l,c) plot imagedata i, with styles p,l,c is for
%  points, lines & conics, respectively.

numberpoints=0;
if nargin>1 & strncmp(p,'numbered',8),
  numberpoints=1;
  p=p(9:end);
end
if nargin<=1 | isempty(p),
  p = '*';
end
if nargin<=2 | isempty(l),
  l = '-';
end
if nargin<=3 | isempty(c),
  c = '.';
end

holdmode=ishold;

if ~isempty(i.im),
  %clf;
  if issparse(i.im)
    axis([0.5,size(i.im,2)+0.5,0.5,size(i.im,1)+0.5]);
    axis('ij');
  else
    image(i.im);
  end
  noimage=0;
elseif ~isempty(i.filename),
  %clf;
  if length(i.filename)>4 & strcmp(i.filename(end-3:end),'.pgm'),
    image(readpgm(i.filename));
  else
    image(imread(i.filename));
  end
  noimage=0;
else
  if holdmode == 0,
    %clf;
    axis('ij');
  end
  noimage=1;
end
hold on;
if noimage == 0,
  colormap(gray(256));
end
ax = axis;

if size(i.points,2)>0,
  plot( (i.points(1,:)./i.points(3,:))',...
	(i.points(2,:)./i.points(3,:))', p);
  if numberpoints,
    pts=[i.points(1,:)./i.points(3,:);...
	 i.points(2,:)./i.points(3,:)];
    ax=axis;dev=0.003*ax(2);
    for pp=1:size(pts,2);
      h=text(double(pts(1,pp)-dev),double(pts(2,pp)-dev),num2str(pp));
      set(h,'Fontsize',11);set(h,'Color','g');
    end
  end

end

plotconic(i.conics,c);
plotline(i.lines,l);

if holdmode == 0
  hold off;
end






