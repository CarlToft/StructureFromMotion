function imseq=project(s,m,v);
% PROJECT imseq=project(s,m,v) projects structure s with motion m
%        in view v to an image sequence
% INPUTS:
%  s - structure
%  m - motion
%  v - views, optional. If not specified all views are assumed
% OUTPUT:
%  imseq - cell of imagedata objects. If only one view, imseq is
%           an imagedata object


if nargin<3,
  v=1:size(m);
end

imseq=cell(1,length(v));
cnt=1;

for view=v;
  P=getcameras(m,view);
  % project points
  U=getpoints(s);
  if ~isempty(U);
   u=P*U;
   ii=find(abs(u(3,:))>1e-5);
   u(1,ii)=u(1,ii)./u(3,ii);
   u(2,ii)=u(2,ii)./u(3,ii);
   u(3,ii)=1;
  else
   u=[];
  end

  % project lines
  L=getlines(s);
  if ~isempty(L);
   u1= P*L(1:4,:);
   u2= P*L(5:8,:);

   %store end points
   l=[u1;u2];
%   l=zeros(size(u1));
%   for i=1:size(u1,2);
%    if ~isnan(u1(:,i))
%      li=null([u1(:,i)';u2(:,i)']);
%    else
%      li=NaN*ones(3,1);
%    end
%    l(:,i)=li;
%   end
   
  else
   l=[];
  end


  %project quadrics
  Q=getquadrics(s);
  if ~isempty(Q);
   Pconic = coniccamera(P);
   c=Pconic*Q;
  else
   c=[];
  end

  imseq{cnt}=imagedata([],u,l,c);
  cnt=cnt+1;
end

if length(v)==1,
  imseq=imseq{1};
end

  
