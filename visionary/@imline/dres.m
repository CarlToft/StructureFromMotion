function [r,drdP,drdU]=dres(im,proj,ob);
% IMLINE/DRES function [r,drdP,drdU]=dres(im,proj,ob);
% calculates weighted residuals and their derivatives with
% respect to changes in structure and motion

[U,dU]=udu(ob);

[P,dP]=pdp(proj);

upl = P*U;


dupldU=zeros(size(upl,1),size(upl,2),sizedx(ob));
for i=1:sizedx(ob);
  dupldU(:,:,i)=P*squeeze(dU(:,:,i));
end;
dupldP=zeros(size(upl,1),size(upl,2),sizedx(proj));
for i=1:sizedx(proj);
  dupldP(:,:,i)=squeeze(dP(:,:,i))*U;
end;


up1=upl(:,1);
up2=upl(:,2);
dup1dU=squeeze(dupldU(:,1,:));
dup2dU=squeeze(dupldU(:,2,:));
dup1dP=squeeze(dupldP(:,1,:));
dup2dP=squeeze(dupldP(:,2,:));



um=im.u;

n=[0,0,1]';
un1=up1/(up1'*n);
un2=up2/(up2'*n);


dun1dup1 = (eye(size(up1,1))/(up1'*n) - up1*n'/(up1'*n)^2);
dun2dup2 = (eye(size(up2,1))/(up2'*n) - up2*n'/(up2'*n)^2);
dun1dP = dun1dup1*dup1dP;
dun1dU = dun1dup1*dup1dU;
dun2dP = dun2dup2*dup2dP;
dun2dU = dun2dup2*dup2dU;

r=[un1,un2]'*um/im.stddevs;

drdP=[um'*dun1dP;um'*dun2dP]/im.stddevs;
drdU=[um'*dun1dU;um'*dun2dU]/im.stddevs;




%up = cross(upl(:,1),upl(:,2));
%dupdU=zeros(3,sizedx(ob));
%for i=1:sizedx(ob);
%  dupdU(:,i)=cross(dupldU(:,1,i),upl(:,2))+cross(upl(:,1),dupldU(:,2,i));
%end;
%dupdP=zeros(3,sizedx(proj));
%for i=1:sizedx(proj);
%  dupdP(:,i)=cross(dupldP(:,1,i),upl(:,2))+cross(upl(:,1),dupldP(:,2,i));
%end;

%n=im.n;
%um=im.u;
%L=im.L;

%un= up/(up'*n);
%dundup = (eye(size(up,1))/(up'*n) - up*n'/(up'*n)^2);
%dundP = dundup*dupdP;
%dundU = dundup*dupdU;

%r = L*(um-un);
%drdP = -L*dundP;
%drdU = -L*dundU;


