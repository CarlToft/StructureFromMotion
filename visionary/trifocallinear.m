function T=trifocallinear(imdata,index,points);
% function T=trifocallinear(imdata,index,points);
% Determines the trifofal tensor linearly using point correspondences
%
% Input:
%   imdata - cell of IMAGEDATA
%   index - indices of images to use
%   points - indices of points to use
% Output:
%   T - trifocal tensor

if nargin<2,
  index = [ 1 2 3];
end;

if nargin<3,
  points=getcommonfeatures(imdata(index));
end;

imT1=getnormtrans(imdata{index(1)});
imT2=getnormtrans(imdata{index(2)});
imT3=getnormtrans(imdata{index(3)});

x1=imT1*getpoints(imdata{index(1)},points);
x2=imT2*getpoints(imdata{index(2)},points);
x3=imT3*getpoints(imdata{index(3)},points);


%kb jb kp jp k j
J=[1 1 3 3 2 2 1
   1 1 3 2 2 3 -1
   1 2 3 3 2 1 -1
   1 2 3 1 2 3 1
   1 3 3 1 2 2 -1
   1 3 3 2 2 1 1
   
   2 1 1 3 3 2 1
   2 1 1 2 3 3 -1
   2 2 1 3 3 1 -1
   2 2 1 1 3 3 1
   2 3 1 1 3 2 -1
   2 3 1 2 3 1 1
      
   3 1 2 3 1 2 1
   3 1 2 2 1 3 -1
   3 2 2 3 1 1 -1
   3 2 2 1 1 3 1
   3 3 2 1 1 2 -1
   3 3 2 2 1 1 1
   
   
   1 1 2 3 3 2 -1
   1 1 2 2 3 3 1
   1 2 2 3 3 1 1
   1 2 2 1 3 3 -1
   1 3 2 1 3 2 1
   1 3 2 2 3 1 -1
   
   2 1 3 3 1 2 -1
   2 1 3 2 1 3 1
   2 2 3 3 1 1 1
   2 2 3 1 1 3 -1
   2 3 3 1 1 2 1
   2 3 3 2 1 1 -1
      
   3 1 1 3 2 2 -1
   3 1 1 2 2 3 1
   3 2 1 3 2 1 1
   3 2 1 1 2 3 -1
   3 3 1 1 2 2 1
   3 3 1 2 2 1 -1
 ];

         

M=zeros(9*size(points,2),27);
for n=1:size(points,2)

for m=1:36
kb=J(m,1);jb=J(m,2);kp=J(m,3);jp=J(m,4);k=J(m,5);j=J(m,6);
for i=1:3        
   M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)=M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)+x1(i,n)*x2(jp,n)*x3(kp,n)*J(m,7);
 end
 end 
 end
 
 [u,s,v]=svd(M);
 t=v(:,27);


for i=1:3
  for j=1:3
    for k=1:3

   T(i,j,k)=t(9*(i-1)+3*(j-1)+k);

    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% i ursprungliga bildkoordinater
m=tritop(T);
P1=getcameras(m,1);
P2=getcameras(m,2);
P3=getcameras(m,3);
m=motion({inv(imT1)*P1,inv(imT2)*P2,inv(imT3)*P3});
T=ptotri(m);

 
