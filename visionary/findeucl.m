function e=findeucl(data1,data2,fixscale);
% e=FINDEUCL(data1,data2) - finds euclidean transformation 'e' using points
% "e(data1) = data2"
% data1, data2 : structure or imagedata objects
% See also: changecsystem

if nargin<3,
    fixscale=0;
end

  x = pflat(getpoints(data1));
  y = pflat(getpoints(data2));
  
  m = size(x,1)-1;
  x = x(1:m,:);
  y = y(1:m,:);

  n = size(x,2);
  X1= y; X2 = x;
  
  m1 = mean(X1');  
  m2 = mean(X2');  
  
  X1 = X1-m1'*ones(1,n); %determining centroid
  X2 = X2-m2'*ones(1,n);
if fixscale,
    s1=1;
    s2=1;
else
  s1 = sum(sum(X1.^2).^.5); %determining scale
  s2 = sum(sum(X2.^2).^.5); 
end

  X1 = X1/s1;
  X2 = X2/s2;
  

  M = zeros(size(X1,1)); %determining rotation
  for i =1:n
    M = M+X1(:,i)*X2(:,i)';
  end
  
  [u,d,v] = svd(M);
  
    

  R = u*v'; 

  if det(R)<0, 
    eye(m); ans(m,m) = -1;
    R = u*ans*v'; 
  end
    e = [s1/s2*R, -s1/s2*R*(m2')+m1'; zeros(1,m) 1];
