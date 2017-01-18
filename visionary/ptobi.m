function B=ptobi(m,index)
% B=ptobi(m,index)
% Calculates the bilinear tensor from the two camera matrices in MOTION m
%   index: specifies camera matrices in m (optional)
%   Default: index = [1 2]

if nargin<2,
  index = [ 1 2];
end;

P1=getcameras(m,index(1));
P2=getcameras(m,index(2));

for i=1:3
  for j=1:3
    B(i,j)=det([P1(1+mod(i,3),:) ; P1(1+mod(i+1,3),:) ; P2(1+mod(j,3),:) ; P2(1+mod(j+1,3),:)]);
  end
end
