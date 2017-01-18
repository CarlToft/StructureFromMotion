function plotline(lines,style)
%PLOTLINE plot(lines,style) plot lines
%  lines - 3xn (8xn) matrix where columns are homog. lines
%  style - style parameters for plot
%  The current axis of figure are used as limitations of the lines
%  2D lines are represented by homogeneous coordinates
%  (or by two homogeneous coordinates for end points).
%  3D lines are represented by two homogeneous 3D points.

if nargin<2,
 style = '-';
end

numberlines=0;
if strncmp(style,'numbered',8),
  numberlines=1;
  style=style(9:end);
  if length(style)==0,
   style = '-';
  end
end


holdmode=ishold;

if size(lines,1)==3,
 ax = axis;
 for t = 1:size(lines,2)
  tmp = lines(:,t);
  if abs(tmp(1))<abs(tmp(2)),
    x1=ax(1);
    y1=(-tmp(3)-tmp(1)*x1)/tmp(2);
    x2=ax(2);
    y2=(-tmp(3)-tmp(1)*x2)/tmp(2);
  else
    y1=ax(3);
    x1=(-tmp(3)-tmp(2)*y1)/tmp(1);
    y2=ax(4);
    x2=(-tmp(3)-tmp(2)*y2)/tmp(1);
  end
  plot ([x1;x2],[y1;y2],style);
  hold on;
 end
elseif size(lines,1)==6, %2d line segment
 ax = axis;
 for t = 1:size(lines,2),
  tmp = lines(:,t);
  xx1=pflat(tmp(1:3));
  xx2=pflat(tmp(4:6));
  plot ([xx1(1);xx2(1)],[xx1(2);xx2(2)],style);
  hold on;
  
  if numberlines,
    ax=axis;dev=0.003*ax(2);
    h=text((xx1(1)+xx2(1))/2-dev,(xx1(2)+xx2(2))/2-dev,num2str(t));
    set(h,'Fontsize',11);set(h,'Color','y');
  end
  
  
 end
elseif size(lines,1)==8, %3D line segment
 if isempty(get(gcf,'CurrentAxes'))
   ax = [0 1 0 1 0 1];
 else
   ax = axis;
   if size(ax,2)==4;
     ax =[ax ax(1:2)];
   end
 end
 for t = 1:size(lines,2)
  if ~isnan(lines(1,t))
    u1 = pflat(lines(1:4,t)); u1 = u1(1:3);
    u2 = pflat(lines(5:8,t)); u2 = u2(1:3);

    p=[u1';u2'];
%    d = normc(u1-u2);
%    v = abs(d);
%    i = find(v==max(v));
%    t1 = (ax(2*i-1)-u1(i))/d(i);
%    t2 = (ax(2*i)  -u1(i))/d(i);
%    p = [u1' + t1*d'; u1' + t2*d'];

    plot3(p(:,1),p(:,2),p(:,3),style);
    hold on;
  end
 end
end

if holdmode == 0
  hold off;
end




