function T=ptotri(m,index)
% T=ptotri(m,index)
% Calculates the trilinear tensor from the camera matrices in MOTION m
%   index: specifies camera matrices in m (optional)
%   Default: index = [1 2 3]

if nargin<2,
  index = [ 1 2 3];
end;

P1=getcameras(m,index(1));P1=P1/norm(P1);
P2=getcameras(m,index(2));P2=P2/norm(P2);
P3=getcameras(m,index(3));P3=P3/norm(P3);


for i=1:3
  for j=1:3
    for k=1:3
      T(i,j,k)=det([P1(1+mod(i,3),:) ; P1(1+mod(i+1,3),:) ; P2(j,:) ; P3(k,:)]);
    end
  end
end
