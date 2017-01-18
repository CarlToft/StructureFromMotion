function pts=getpointsconic(u, conicsampling)
% WORKS ONLY FOR ELLIPSE
  A  = [u(1), u(2)/2; u(2)/2, u(3)];
  bb = [u(4); u(5)];
  c  = u(6);
  [Q D] = eig(A);
  det   = D(1,1)*D(2,2);
  if (det == 0),
%error('not implemented');
pts=[];
  else
    bs    = Q'*bb;
    alpha = atan2(Q(2,1), Q(1,1));
    zs    = -(2*D)\bs;  
    z     = Q*zs;
    h     = -bs'*zs/2-c;
    a     = h/D(1,1);
    b     = h/D(2,2);
    if ((a > 0) & (b > 0)),
      pts=getellipse (z, sqrt(a), sqrt(b), alpha, conicsampling);
    else
pts=[];
%error('not implemented');
    end
  end % if

%end % drawconic

function pts=getellipse (z, a, b, alpha, conicsampling)
%
%       z, a, b, alpha: parameters of ellipse

  s = sin(alpha); c = cos(alpha);
  Q =[c -s; s c];
  theta = [0:2*pi/conicsampling:2*pi];
  pts = diag(z)*ones(2,length(theta)) + ...
      Q*[a*cos(theta); b*sin(theta)];

  pts(:,size(pts,2))=[];
%end drawellipse





