function plotconic(conics,style)
%PLOTCONIC plot conics - 6xn matrix where columns are dual conics

if nargin<2,
 style = '-';
end

for t = 1:size(conics,2)
  v = pflat(conics(:,t));
  if ~isnan(v(1));
    c = inv(v2m(v));
    drawconic([c(1,1),c(1,2)*2,c(2,2),c(1,3)*2,c(2,3)*2,c(3,3)],style);
  end
end

function drawconic (u,style)
%DRAWCONIC
%
%       drawconic (u,style)
%       draw solutions of x'Ax + bb'x + c = 0.
%
  A  = [u(1), u(2)/2; u(2)/2, u(3)];
  bb = [u(4); u(5)];
  c  = u(6);

  [Q D] = eig(A);
  det   = D(1,1)*D(2,2);
  if (det == 0),
    if ((D(1,1) == 0) & (D(2,2) == 0)),
      z     = -c*(bb' / norm(bb, 2));
      alpha = atan2(-bb(1), bb(2));
      drawline (z, alpha,style);
    else
      if (D(1,1) == 0),
        D = D*[0 1; 1 0];
        Q = Q*[0 1; 1 0];
      end
      alpha = atan2(Q(2,1), Q(1,1));
      z(1)  = bb(1)/2/D(1,1);
      a     = D(1,1)/bb(2);
      z(2)  = (z(1)^2 - c)*a;
      drawparabola (z, a, alpha, style);
    end
  else
    bs    = Q'*bb;
    alpha = atan2(Q(2,1), Q(1,1));
    zs    = -(2*D)\bs;  
    z     = Q*zs;
    h     = -bs'*zs/2-c;
    a     = h/D(1,1);
    b     = h/D(2,2);
    if ((a > 0) & (b > 0)),
      drawellipse (z, sqrt(a), sqrt(b), alpha, style);
    else
      if (a < 0),
        tmp = a;
        a = b;
        b = -a;
        alpha = alpha + pi/2;
      else
        b = -b;
      end
      if b>=0,
        drawhyperbola (z, sqrt(a), sqrt(b), alpha, style);
      end
    end
  end % if

%end % drawconic

function drawline (C, alpha, style)
%DRAWLINE
%
%       drawline(C, alpha, style)
%       draws line through C.

  s = sin(alpha); c = cos(alpha);
  Q =[c -s; s c];
  theta = [-10:0.2:10];
  u = diag(C)*ones(2,length(theta)) + ...
       Q*[sin(theta); cos(theta)];
  plot(u(1,:),u(2,:),style);

%end % drawline


function drawparabola (z, a, alpha, style)
%DRAWPARABOLA
%
%       drawparabola (z, a, alpha, style)
%       draw parabola for a*(x - z(1))^2 + (y - z(2)) = 0
%       rotated by alpha
%
  s = sin(alpha); c = cos(alpha);
  Q = [c -s; s c];
  theta = [-10:0.2:10];
  u = diag(z)*ones(2,length(theta)) + ...
      Q*[theta; a*theta.*theta];
  plot(u(1,:),u(2,:), style);

%end % drawparabola

function drawellipse (z, a, b, alpha, style)
%DRAWELLIPSE    Draw ellipse
%
%       drawellipse (z, a, b, alpha, style)
%       draws ellipse into current figure.
%
%       z, a, b, alpha: parameters of ellipse
%       style: pattern to be used

  s = sin(alpha); c = cos(alpha);
  Q =[c -s; s c];
  theta = [0:0.01:2*pi];
  u = diag(z)*ones(2,length(theta)) + ...
      Q*[a*cos(theta); b*sin(theta)];
  plot(u(1,:),u(2,:), style);
%end drawellipse

function drawhyperbola(C, a, b, alpha, style)
%DRAWHYPERBOLA
%
%       drawhyperbola(C, a, b, alpha, style)
%       draw hyperbola with center C  semiaxis a  and b
%       alpha angle between a and x-axe
%
  s = sin(alpha); c = cos(alpha);
  Q = [c -s; s c];
  fac = (a + b + 2)*log(a + b + 2);
  theta = [-fac:fac / 50:fac];
  u = diag(C)*ones(2,length(theta)) + ...
      Q*[a*cosh(theta); b*sinh(theta)];
  plot(u(1,:),u(2,:),style);
  u = diag(C)*ones(2,length(theta)) + ...
      Q*[-a*cosh(theta); b*sinh(theta)];
  plot(u(1,:),u(2,:),style);
%end % drawhyperbola







