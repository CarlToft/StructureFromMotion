function [bestmot,ind,bestd]=ransac2views(imdata,iteration,threshold)
% [mot,ind]=ransac2views(imdata,iteration,threshold)
% Calculate the motion of two views using ransac
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
  iteration=500;
end;
if nargin<3
   threshold=1;
end

xx1=getpoints(imdata{1});
xx2=getpoints(imdata{2});



n1=size(xx1,2);
consistent=0;

M=zeros(7,9);


for i=1:iteration,

  Per=randperm(n1);
  
  for jj=1:7,
      index=Per(jj);
      x1=xx1(1,index);
      x2=xx1(2,index);
      x3=xx2(1,index);
      x4=xx2(2,index);
      M(jj,:)=[x1*x3,x2*x3,x3,x1*x4,x2*x4,x4,x1,x2,1];
  end
  [u,s,v]=svd(M);
  f1=v(1,8);f2=v(2,8);f3=v(3,8);f4=v(4,8);f5=v(5,8);f6=v(6,8);f7=v(7,8);f8=v(8,8);f9=v(9,8);
  g1=v(1,9);g2=v(2,9);g3=v(3,9);g4=v(4,9);g5=v(5,9);g6=v(6,9);g7=v(7,9);g8=v(8,9);g9=v(9,9);
  F=[f1,f2,f3;f4,f5,f6;f7,f8,f9];
  G=[g1,g2,g3;g4,g5,g6;g7,g8,g9];
  
  poly=[g1*g5*g9-g1*g6*g8-g4*g2*g9+g4*g3*g8+g7*g2*g6-g7*g3*g5,g1*f5*g9+g1*g5*f9-g1*f6*g8-g1*g6*f8+f1*g5*g9-f1*g6*g8-g4*f2*g9-g4*g2*f9+g4*f3*g8+g4*g3*f8-f4*g2*g9+f4*g3*g8+g7*f2*g6+g7*g2*f6-g7*f3*g5-g7*g3*f5+f7*g2*g6-f7*g3*g5,...
          g1*f5*f9-g1*f6*f8+f1*f5*g9+f1*g5*f9-f1*f6*g8-f1*g6*f8-g4*f2*f9+g4*f3*f8-f4*f2*g9-f4*g2*f9+f4*f3*g8+f4*g3*f8+g7*f2*f6-g7*f3*f5+f7*f2*g6+f7*g2*f6-f7*f3*g5-f7*g3*f5,...
          f4*f3*f8-f7*f3*f5+f7*f2*f6-f4*f2*f9-f1*f6*f8+f1*f5*f9];
  ll=roots(poly);
  
  ll=ll(find(abs(imag(ll))==0));

  for jj=1:length(ll);
      
      H=F+ll(jj)*G;
      
      tmp=log(svd(H/norm(H))+eps);
      if tmp(2)>-5, %rank one matrix not allowed
          
        %check feasibility of other points
        % (HARTLEY/STURM triangulation)
        H=H'/norm(H);
        e1=pflat(cross(H(:,1),H(:,2)));
        e2=pflat(cross(H(1,:)',H(2,:)'));

        T1=eye(3);T1(:,3)=e1;
        T2=eye(3);T2(:,3)=e2;

        tmp=T1'*H*T2;
        [u,s,v]=svd(tmp(1:2,1:2));
        T1(1:2,1:2)=T1(1:2,1:2)*u';
        T2(1:2,1:2)=T2(1:2,1:2)*v;
        Ftmp=T1'*H*T2;
        lambda=Ftmp(1,1)/Ftmp(2,2);

        xx1t=inv(T1)*xx1;
        xx2t=inv(T2)*xx2;

        res=zeros(size(xx1t,2),1);
        for qq=1:size(xx1t,2),
            mx1=xx1t(1,qq);
            my1=xx1t(2,qq);
            mx2=xx2t(1,qq);
            my2=xx2t(2,qq);
    

            poly=[mx1*my1-mx2*lambda*my2,...
    -lambda^2*my2^2+mx2^2*lambda^2-mx1^2+my1^2,...
    -mx1*my1+2*mx1*my1*lambda^2+mx2*lambda^3*my2-2*mx2*lambda*my2,...
    -2*lambda^2*my2^2-2*mx1^2*lambda^2+2*mx2^2*lambda^2+2*my1^2*lambda^2,...
    mx1*my1*lambda^4-2*mx1*my1*lambda^2+2*mx2*lambda^3*my2-mx2*lambda*my2,...
    -mx1^2*lambda^4-lambda^2*my2^2+mx2^2*lambda^2+my1^2*lambda^4,...
    -mx1*my1*lambda^4+mx2*lambda^3*my2];

            r=roots(poly);
            r=r(find(imag(r)==0));

            [dist,indqq]=min((r*mx1+my1).^2./(r.^2+1)+(lambda*mx2-r*my2).^2./(r.^2+lambda^2));
            res(qq)=dist;
        end

        
      
        index=find(res<threshold^2);
      
        if length(index)>consistent,
          consistent=length(index);
          ind=index;
          bestd=res;
          bestF=H;
        end
    end %if rank one
  end
end

bestmot=bitop(bestF);


