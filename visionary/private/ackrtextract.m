function  [K,R,t]=ackrtextract(m);
nbrcam=length(m);
K=zeros(nbrcam,5);
R=cell(nbrcam,1);
t=cell(nbrcam,1);

for i=1:nbrcam
 P =pdp(m{i});
 [KK,Rt]=rq(P);KK=KK/KK(3,3);
 K(i,1)=KK(1,1);
 K(i,2)=KK(2,2)/KK(1,1);
 K(i,3)=KK(1,2)/KK(1,1);
 K(i,4)=KK(1,3);
 K(i,5)=KK(2,3);
 R{i}=Rt(:,1:3);
 t{i}=-R{i}'*Rt(:,4);
end

 
