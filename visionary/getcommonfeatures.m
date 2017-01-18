function [p,l,c]=getcommonfeatures(imseq);
% [p,l,c]=getcommonfeatures(imseq) returns indexes of point, line
%   and conic features present in whole sequence

p0=ones(1,size(imseq{1},1));
l0=ones(1,size(imseq{1},2));
c0=ones(1,size(imseq{1},3));

for i=1:length(imseq);

 pi=getpoints(imseq{i});
 if size(pi,2)>0
   pi=~isnan(pi(1,:));
   p0 = p0 & pi;
 else
   p0=0;
 end

 li=getlines(imseq{i});
 if size(li,2)>0
   li=~isnan(li(1,:));
   l0 = l0 & li;
 else
   l0=0;
 end

 ci=getconics(imseq{i});
 if size(ci,2)>0
   ci=~isnan(ci(1,:));
   c0 = c0 & ci;
 else
   c0=0;
 end
end

p=find(p0);
l=find(l0);
c=find(c0);
