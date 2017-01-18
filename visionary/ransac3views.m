function [bestmot,ind,bestd]=ransac3views(imdata,iteration,threshold)
% [mot,ind]=ransac3views(imdata,iteration,threshold)
% Calculate the motion of three views using ransac
% Input:
%   imdata - cell of IMAGEDATA
%   iteration - iterations in RANSAC (optional)
%   threshold - for inliers in pixels (optional)
% Output:
%   mot - motion
%   ind - indices of used points
%   bestd - reprojection errors
% Default: iteration=100, threshold=1;

if nargin<2,
  iteration=100;
end;
if nargin<3
   threshold=1;
end

swarn=warning;
warning('off');

imT1=getnormtrans(imdata{1});
imT2=getnormtrans(imdata{2});
imT3=getnormtrans(imdata{3});
invT1=inv(imT1);
invT2=inv(imT2);
invT3=inv(imT3);


x1orig=getpoints(imdata{1});
x2orig=getpoints(imdata{2});
x3orig=getpoints(imdata{3});


x1all=psphere(imT1*x1orig);
x2all=psphere(imT2*x2orig);
x3all=psphere(imT3*x3orig);

imseqnew={imagedata([],x1all),imagedata([],x2all),imagedata([],x3all)};


n1=size(x1all,2);
consistent=0;


M=zeros(9,7);


for i=1:iteration,

  Per=randperm(n1);
  Per=Per(1:6);

  [scell,mcell]=sm6points(imseqnew,Per);

  for j=1:length(scell);

    P1=getcameras(mcell{j},1);
    P2=getcameras(mcell{j},2);
    P3=getcameras(mcell{j},3);
    
    tmp1=svd(P1);
    tmp2=svd(P2);
    tmp3=svd(P3);
    if tmp1(3)>1e-8 & tmp2(3)>1e-8 & tmp3(3)>1e-8,
    

     M(1:3,1:4)=P1/norm(P1);
     M(4:6,1:4)=P2/norm(P2);
     M(7:9,1:4)=P3/norm(P3);

     P1orig=invT1*P1;
     P2orig=invT2*P2;
     P3orig=invT3*P3;


     p1a=P1orig(1,:);
     p2a=P1orig(2,:);
     p3a=P1orig(3,:);
     p1b=P2orig(1,:);
     p2b=P2orig(2,:);
     p3b=P2orig(3,:);
     p1c=P3orig(1,:);
     p2c=P3orig(2,:);
     p3c=P3orig(3,:);
     d=zeros(1,n1);

     for qq=1:n1;
      if all(Per-qq)
	M(1:3,5)=x1all(:,qq);
	M(4:6,6)=x2all(:,qq);
	M(7:9,7)=x3all(:,qq);
	[u,s,v]=svd(M);
	X=v(1:4,7)/v(4,7);

        for rr=1:2, %two bundlesteps
	  p3aX=p3a*X;p3bX=p3b*X;p3cX=p3c*X;
          p1aX=p1a*X/p3aX;p2aX=p2a*X/p3aX;
          p1bX=p1b*X/p3bX;p2bX=p2b*X/p3bX;
          p1cX=p1c*X/p3cX;p2cX=p2c*X/p3cX;
	  f=[p1aX-x1orig(1,qq);p2aX-x1orig(2,qq);
	     p1bX-x2orig(1,qq);p2bX-x2orig(2,qq);
	     p1cX-x3orig(1,qq);p2cX-x3orig(2,qq)];
          fgrad=[p1a(1:3)/p3aX-p1aX/p3aX*p3a(1:3);
		 p2a(1:3)/p3aX-p2aX/p3aX*p3a(1:3);
		 p1b(1:3)/p3bX-p1bX/p3bX*p3b(1:3);
		 p2b(1:3)/p3bX-p2bX/p3bX*p3b(1:3);
		 p1c(1:3)/p3cX-p1cX/p3cX*p3c(1:3);
		 p2c(1:3)/p3cX-p2cX/p3cX*p3c(1:3)];
	  X(1:3)=X(1:3)-inv(fgrad'*fgrad+0.0001*eye(3))*fgrad'*f;
        end
        p3aX=p3a*X;p3bX=p3b*X;p3cX=p3c*X;
        d(qq)=(p1a*X/p3aX-x1orig(1,qq))^2+(p2a*X/p3aX-x1orig(2,qq))^2+...
              (p1b*X/p3bX-x2orig(1,qq))^2+(p2b*X/p3bX-x2orig(2,qq))^2+...
              (p1c*X/p3cX-x3orig(1,qq))^2+(p2c*X/p3cX-x3orig(2,qq))^2;
      end
    end
    newconsistent=sum(d<15*threshold);
    if newconsistent>consistent;
      consistent=newconsistent;
      ind=find(d<15*threshold);
      bestd=d;
      bestmot=motion({P1orig/norm(P1orig),P2orig/norm(P2orig),P3orig/norm(P3orig)});
    end
    
   end %if singular values of camera matrices...
  end %loop over the "3" solutions

end



warning(swarn);
