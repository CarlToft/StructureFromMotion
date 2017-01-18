function d=taylorerror(imdata,index,T,ptsindex);

if nargin<4,
  ptsindex=1:size(imdata{1},1);
end;

 
imT1=getnormtrans(imdata{index(1)});
imT2=getnormtrans(imdata{index(2)});
imT3=getnormtrans(imdata{index(3)});

x1=imT1*getpoints(imdata{index(1)},ptsindex);
x2=imT2*getpoints(imdata{index(2)},ptsindex);
x3=imT3*getpoints(imdata{index(3)},ptsindex);

m=tritop(T);
P1=getcameras(m,1);
P2=getcameras(m,2);
P3=getcameras(m,3);
m=motion({imT1*P1,imT2*P2,imT3*P3});
T=ptotri(m);



for i=1:3
  for j=1:3
    for k=1:3

   t(9*(i-1)+3*(j-1)+k)=T(i,j,k);

    end
  end
end
t=t';

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
 
n1=size(x1,2);
n2=size(x2,2);
n3=size(x3,2);


M=zeros(9*n1,27);
for n=1:n1
for m=1:36
kb=J(m,1);jb=J(m,2);kp=J(m,3);jp=J(m,4);k=J(m,5);j=J(m,6);
for i=1:3        
   M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)=M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)+x1(i,n)*x2(jp,n)*x3(kp,n)*J(m,7);
 end
 end 
 end

R=M*t;
 
for i=1:n1

r=R(9*(i-1)+1:9*i);
jJ=trijacob(t,[x1(:,i),x2(:,i),x3(:,i)]);
d(i)=r'*pinv(jJ*jJ')*r;

end
dd=abs(eig(imT1));
scale=mean(dd(1:2))/dd(3);
d=sqrt(d/scale/6)*10;




%%%%%%%%%%%%%%%%%%%%%%%%%%


function jJ=trijacob(t,points)

mm1=points(:,1);
mm2=points(:,2);
mm3=points(:,3);
m1=mm1;
m2=mm2;
m3=mm3;
n=1;

 
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

 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
m1(:,1)=[1 0 0]';

M=zeros(9,27);
for m=1:36
kb=J(m,1);jb=J(m,2);kp=J(m,3);jp=J(m,4);k=J(m,5);j=J(m,6);
for i=1:3       
   M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)=M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)+m1(i,n)*m2(jp,n)*m3(kp,n)*J(m,7);
 end
 end

 
jJ(:,1)=M(1:9,:)*t;
m1=mm1;
%%%%%%%%%%%%%%%%%%%%%%

m1(:,1)=[0 1 0]';

M=zeros(9,27);
for m=1:36
kb=J(m,1);jb=J(m,2);kp=J(m,3);jp=J(m,4);k=J(m,5);j=J(m,6);
for i=1:3       
   M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)=M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)+m1(i,n)*m2(jp,n)*m3(kp,n)*J(m,7);
 end
 end
 
jJ(:,2)=M(1:9,:)*t;
m1=mm1;
%%%%%%%%%%%%%%%%%%%%%%

m2(:,1)=[1 0 0]';

M=zeros(9,27);
for m=1:36
kb=J(m,1);jb=J(m,2);kp=J(m,3);jp=J(m,4);k=J(m,5);j=J(m,6);
for i=1:3       
   M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)=M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)+m1(i,n)*m2(jp,n)*m3(kp,n)*J(m,7);
 end
 end
 
jJ(:,3)=M(1:9,:)*t;
m2=mm2;
%%%%%%%%%%%%%%%%%%%%%%

m2(:,1)=[0 1 0]';

M=zeros(9,27);
for m=1:36
kb=J(m,1);jb=J(m,2);kp=J(m,3);jp=J(m,4);k=J(m,5);j=J(m,6);
for i=1:3       
   M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)=M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)+m1(i,n)*m2(jp,n)*m3(kp,n)*J(m,7);
 end
 end
 
jJ(:,4)=M(1:9,:)*t;
m2=mm2;
%%%%%%%%%%%%%%%%%%%%%%

m3(:,1)=[1 0 0]';

M=zeros(9,27);
for m=1:36
kb=J(m,1);jb=J(m,2);kp=J(m,3);jp=J(m,4);k=J(m,5);j=J(m,6);
for i=1:3       
   M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)=M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)+m1(i,n)*m2(jp,n)*m3(kp,n)*J(m,7);
 end
 end
 
jJ(:,5)=M(1:9,:)*t;

m3=mm3;
%%%%%%%%%%%%%%%%%%%%%%

m3(:,1)=[0 1 0]';

M=zeros(9,27);
for m=1:36
kb=J(m,1);jb=J(m,2);kp=J(m,3);jp=J(m,4);k=J(m,5);j=J(m,6);
for i=1:3       
   M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)=M(9*(n-1)+3*(jb-1)+kb,9*(i-1)+3*(j-1)+k)+m1(i,n)*m2(jp,n)*m3(kp,n)*J(m,7);
 end
 end

 
jJ(:,6)=M(1:9,:)*t;

m3=mm3;


 
