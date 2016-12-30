function [str,mot,imseq]=projmotion(filenames,logfile);
% [str,mot]=projmotion(filenames,logfile);
% Calculates projective structure, motion and point correspondences
% from a sequence of closely spaced images
% Input:
%   filenames - Filenames of images
%   logfile - Save log-data during running (optional)
% Output:
%   str - structure
%   mot - motion
%   imseq - cell of IMAGEDATA with points

if nargin<2,
  logfile='';
end


%%%%%%%%%%%% Parameters
ransaciterations=2000;
ransacpoints=6;
threshold=3;
mincornerdistance=30; %pixels
intensitythreshold=3;
patchsize=5;
stdev=1.5;
minpoints=10;
borderignore=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np=2*patchsize+1;
nbrimages=length(filenames);


images=cell(1,nbrimages);
for ii=1:nbrimages,
  if strcmp(filenames{ii}(end-3:end),'.pgm'),
    images{ii}=readpgm(filenames{ii});
  else
    images{ii}=imread(filenames{ii});
    if length(size(images{ii}))==3,
      images{ii}=rgb2gray(images{ii});
    end
  end
end
[szimage1,szimage2]=size(images{1});

%%%%%%%%%%%%%%%%%%%%%
% Determine correspondences in views 1,2,3
%%%%%%%%%%%%%%%%%%%%%

imharris=harris(double(images{2}),{'noimage','threshold=1','mindistance=15'});
%imharris=harris(double(images{2}),{'noimage','threshold=0.7','mindistance=15'});

% For each corner in image 2, find it in images 1 & 3

pos2=getpoints(imharris);
[pos1,res1]=findpatch3(double(images{2}),double(images{1}),pos2);
tmp=find(res1<intensitythreshold);
pos1=pos1(:,tmp);
pos2=pos2(:,tmp);

[pos3,res3]=findpatch3(double(images{2}),double(images{3}),pos2);

tmp=find(res3<intensitythreshold);
pos1=pos1(:,tmp);
pos2=pos2(:,tmp);
pos3=pos3(:,tmp);

imseq={imagedata([],pos1),...
       imagedata([],pos2),...
       imagedata([],pos3)};

%%%%%%%%%%%%%%%%% - end of correspondences 1,2,3
if size(imseq{1},1)<=minpoints,
  disp(['Too few points after views 1-2-3']);
  keyboard;
end

%% Ransac for correct motion
ransacpts=size(imseq{1},1);

%[bestH1,bestH2,indhomo,reshomo]=ransac3homos(imseq,ransaciterations,threshold);

[mot,indransac]=ransac3views(imseq,ransaciterations,threshold);

for i=1:3;
  imseq{i}=imagedata(filenames{i},getpoints(imseq{i},indransac));
end;
%figure(1);plot(imseq{1},'numbered');figure(2);plot(imseq{2},'numbered');figure(3);plot(imseq{3},'numbered');

if size(imseq{1},1)<=minpoints,
  disp(['Too few points after ransac in views 1-2-3']);
  keyboard;
end

str=intsecpoints(imseq,mot);

% bundle for view 1,2,3

[str,mot] = bundleplc(str,mot,imseq,{'iteration=50','outputoff','lambda=1000'});
[str,mot] = bundleplc(str,mot,imseq,{'iteration=50','outputoff','lambda=500'});
%determine which points to keep...
err=zeros(1,size(str,1));
for i=1:3;
  imreproj=project(str,mot,i);
  pts=pflat(getpoints(imseq{i}));
  repts=pflat(getpoints(imreproj));
  err=err+sum((pts-repts).^2);
end
ind=find(err<2*threshold);
for i=1:3;
  imseq{i}=imagedata([],getpoints(imseq{i},ind));
end;
str=structure(getpoints(str,ind));

restri=sum(err(ind));

[str,mot] = bundleplc(str,mot,imseq,{'iteration=50','outputoff','lambda=500'});

%rmspoints(str,mot,imseq);
%reproject(str,mot,imseq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda1=log(6);
lambda2=log(6*ransacpts);
lambda3=2;
sigma=1;
r=6;

%homography
d=2;
k=16;


gric_h=reshomo/sigma^2+(ransacpts-length(indhomo))*lambda3*(r-d)+...
     lambda1*d*ransacpts+lambda2*k;

%trifocal
d=3;
k=18;

gric_t=restri/sigma^2+(ransacpts-length(indransac))*lambda3*(r-d)+...
     lambda1*d*ransacpts+lambda2*k;

if gric_h<gric_t,
  disp('Homography better...');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find more matches using known motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F12=ptobi(mot,[1,2]);F12=F12/norm(F12);
F23=ptobi(mot,[2,3]);F23=F23/norm(F23);
F13=ptobi(mot,[1,3]);F13=F13/norm(F13);

%examine "unused" harris-corners

imtmp=imagedata;
harrispts=getpoints(imharris);
pts=getpoints(imseq{2});
for i=1:size(harrispts,2);
  dist=(harrispts(1,i)-pts(1,:)).^2+(harrispts(2,i)-pts(2,:)).^2;
  if min(dist)>mincornerdistance^2,
    imtmp=addpoints(imtmp,harrispts(:,i));
  end
end

newimseq=findmorematches(images(1:3),imtmp,F12,F23,F13,threshold,4*intensitythreshold/3);

newstr=intsecpoints(newimseq,mot);
newstr=bundleplc(newstr,mot,newimseq,{'lambda=1','iteration=50','outputoff','structure'});


%calculate reprojection error
err=zeros(1,size(newstr,1));
for i=1:3;
  imreproj=project(newstr,mot,i);
  pts=pflat(getpoints(newimseq{i}));
  repts=pflat(getpoints(imreproj));
  err=err+sum((pts-repts).^2);
end
ind=find(err<15*threshold);


for i=1:3;
  imseq{i}=imseq{i}+imagedata([],getpoints(newimseq{i},ind));
end;
str=str+structure(getpoints(newstr,ind));

[str,mot] = bundleplc(str,mot,imseq,{'iteration=50','outputoff','lambda=500'});

%FINALLY determine which points to keep...
err=zeros(1,size(str,1));
for i=1:3;
  imreproj=project(str,mot,i);
  pts=pflat(getpoints(imseq{i}));
  repts=pflat(getpoints(imreproj));
  err=err+sum((pts-repts).^2);
end
ind=find(err<2*threshold);
for i=1:3;
  imseq{i}=imagedata([],getpoints(imseq{i},ind));
end;
str=structure(getpoints(str,ind));

[str,mot] = bundleplc(str,mot,imseq,{'iteration=50','outputoff','lambda=500'});


%%%%%%%%%
nbrpoints=size(imseq{1},1);
%originalpatch=2*ones(1,nbrpoints);

if length(logfile)>0,
  eval(['save ',logfile,num2str(3),' str mot imseq nbrpoints']);
end;

%save dummy

%for q=1:3;figure(q);clf;plot(imagedata(double(images{q}),getpoints(imseq{q})));zoom on;end


for imnr=4:nbrimages;

  disp(['Image ',num2str(imnr),' Points: ',num2str(nbrpoints)]);

  prevpts=getpoints(imseq{imnr-1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%follow old paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%


  [currpts,reslist]=findpatch3(double(images{imnr-1}),double(images{imnr}),prevpts,[],130);
  tmp=find(reslist>intensitythreshold);
  currpts(:,tmp)=NaN;




%  im1rho=conv2(double(images{imnr}),ones(np),'valid');
%  im1rho(find(im1rho==0))=1;
%  im1rho=1./im1rho;
%  im1squared=conv2(double(images{imnr}).^2,ones(np),'valid');

%  for ii=1:size(prevpts,2);
%    if ~isnan(prevpts(1,ii)),

%     index=originalpatch(ii);
%      index=imnr-1;
%      inpos=getpoints(imseq{index},ii);
%      predpos=prevpts(:,ii);
%      [newpos,res]=findpatch3(double(images{index}),double(images{imnr}),inpos,predpos,im1rho,im1squared);
%      if res<intensitythreshold,
%	currpts(:,ii)=newpos;
%      end
%    end
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find new corners
%%%%%%%%%%%%%%%%%%%%%%%%%%%

  imharris=harris(double(images{imnr-1}),{'noimage','threshold=1','mindistance=30'});
%  imharris=harris(double(images{imnr-1}),{'noimage','threshold=0.8','mindistance=15'});

  F=ptobi(mot,[imnr-1,imnr-2]);F=F/norm(F);

  plist2=[];
  for ii=1:size(imharris,1);
    pt=getpoints(imharris,ii);
    tmp=prevpts-pt*ones(1,size(prevpts,2));
    dist=tmp(1,:).^2+tmp(2,:).^2;
    [mindist,minindex]=min(dist);
    if mindist>mincornerdistance^2 | isnan(currpts(1,minindex)),
      plist2=[plist2,pt];
    end
  end
  [plist1,reslist1]=findpatchepi(double(images{imnr-1}),double(images{imnr-2}),plist2,F);
  pind1=find(reslist1<intensitythreshold);
  plist1=plist1(:,pind1);
  plist2=plist2(:,pind1);

  [plist3,reslist3]=findpatch3(double(images{imnr-1}),double(images{imnr}),plist2,[],130);
  pind3=find(reslist3<intensitythreshold);

  plist1=plist1(:,pind3);
  plist2=plist2(:,pind3);
  plist3=plist3(:,pind3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of find new corners
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
% Ransac to get new camera matrix
%%%%%%%%%%%%%%%%%%%%

  %create temporary image sequence for imnr-2:imnr

  goodindex=find(~isnan(currpts(1,:)));

  x1orig=[getpoints(imseq{imnr-2},goodindex),plist1];
  x2orig=[getpoints(imseq{imnr-1},goodindex),plist2];
  x3orig=[currpts(:,goodindex),plist3];

  tmpimseq={imagedata([],x1orig),...
          imagedata([],x2orig),...
          imagedata([],x3orig)};

  [randmot,indransac]=ransac3views(tmpimseq,ransaciterations,threshold);

  % add new points to list
  n0=length(goodindex);

  ind1=find(indransac<=n0);
  ind1=goodindex(indransac(ind1));
  ind2=find(indransac>n0);
  ind2=indransac(ind2)-n0;

  tmp=NaN*ones(size(currpts));
  tmp(:,ind1)=currpts(:,ind1);
  imseq{imnr}=imagedata(filenames{imnr},tmp);

  nbraddpts=length(ind2);
  for qqq=1:imnr,
	if qqq==imnr-2,
	  imseq{qqq}=addpoints(imseq{qqq},plist1(:,ind2));
	elseif qqq==imnr-1,
	  imseq{qqq}=addpoints(imseq{qqq},plist2(:,ind2));
	elseif qqq==imnr,
	  imseq{qqq}=addpoints(imseq{qqq},plist3(:,ind2));
	else
	  imseq{qqq}=addpoints(imseq{qqq},NaN*ones(3,nbraddpts));
	end
  end

  %update structure & motion
	  
% works not so good
%  tmp=[getcameras(randmot,1);getcameras(randmot,2)];
%  tmp2=[getcameras(mot,imnr-2);getcameras(mot,imnr-1)];
%  H=tmp\tmp2;
%  bestP3=getcameras(randmot,3)*H;bestP3=bestP3/norm(bestP3);

% works not so good
%  bestP3=reseccam3(motion(getcameras(mot,imnr-2:imnr-1)),tmpimseq,1:3,indransac);

  tmpstr=structure(getpoints(str,ind1));
  tmpimseq={imagedata([],currpts(:,ind1))};

  bestP3=reseclinear(tmpimseq,tmpstr);

  [dummy,tmpmot]=bundleplc(tmpstr,bestP3,tmpimseq,{'iteration=50','lambda=1e4','outputoff','motion'});


  mot=mot+tmpmot;

  tmpimseq{1}=imagedata(sparse(szimage1,szimage2),x1orig(:,ind2+n0));
  tmpimseq{2}=imagedata(sparse(szimage1,szimage2),x2orig(:,ind2+n0));
  tmpimseq{3}=imagedata(sparse(szimage1,szimage2),x3orig(:,ind2+n0));
  tmpmot=motion(getcameras(mot,imnr-2:imnr));
  tmpstr=intsecpoints(tmpimseq,tmpmot);
  tmpstr=bundleplc(tmpstr,tmpmot,tmpimseq,{'lambda=1','iteration=50','outputoff','structure'});
  str=str+tmpstr;

  % REINITIALIZE due to homography in first triplet...
%  if imnr==4 & gric_h<gric_t,
%    [str,mot]=smshape(imseq);
%    str = bundleplc(str,mot,imseq,{'iteration=40','outputoff','lambda=1e4','structure'});
%    [str,mot] = bundleplc(str,mot,imseq,{'iteration=20','outputoff','lambda=1e3'});
%  end

  %%%%%%%%%%%%%%%%%%%%
  % Bundle on all views
  %%%%%%%%%%%%%%%%%%%%
  [str,mot] = bundleplc(str,mot,imseq,{'iteration=20','outputoff','lambda=1e4'});
  [str,mot] = bundleplc(str,mot,imseq,{'iteration=10','outputoff','lambda=5e3'});


  %Determine which points to keep...
  tmp=pflat(getpoints(project(str,mot,imnr),1:nbrpoints)); 
  index=find(sum((tmp-getpoints(imseq{imnr},1:nbrpoints)).^2)<2*threshold);
  tmp=[NaN*ones(3,nbrpoints),plist3(:,ind2)];tmp(:,index)=currpts(:,index);
  imseq{imnr}=imagedata(filenames{imnr},tmp);
  
  [str,mot] = bundleplc(str,mot,imseq,{'iteration=2','outputoff','lambda=5e3'});

  err=zeros(1,size(str,1));
  nbrappeared=zeros(1,size(str,1));
  for i=1:imnr;
    imreproj=project(str,mot,i);
    pts=pflat(getpoints(imseq{i}));
    repts=pflat(getpoints(imreproj));
    tmp=sum((pts-repts).^2);
    ind=find(~isnan(tmp));
    err(ind)=err(ind)+tmp(ind);
    nbrappeared(ind)=nbrappeared(ind)+1;
  end
  ind=find(err./nbrappeared<2*threshold);
  for i=1:imnr;
    imseq{i}=imagedata(filenames{i},getpoints(imseq{i},ind));
  end;
  str=structure(getpoints(str,ind));

  [str,mot] = bundleplc(str,mot,imseq,{'iteration=5','outputoff','lambda=1e4'});

  if sum(~isnan(sum(getpoints(imseq{imnr}))))<=minpoints,
    disp(['Too few points in image ',num2str(imnr)]);
    keyboard;
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find more matches using known motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  prevpts=getpoints(imseq{imnr-1});
  currpts=getpoints(imseq{imnr});
  lostindex=find(isnan(currpts(1,:))&~isnan(prevpts(1,:)));

  for i=1:length(lostindex);

    pt=getpoints(imseq{imnr-1},lostindex(i));
    destpos=pflat(getcameras(mot,imnr)*getpoints(str,lostindex(i)));

    if destpos(1)>borderignore & destpos(1)<szimage2-borderignore & ...
       destpos(2)>borderignore & destpos(2)<szimage1-borderignore,
%     [destpos2,residual]=affineoptimize(double(images{imnr-1}),double(images{imnr}), ...
%           pt, destpos, patchsize, stdev);

      [destpos2,residual]=findpatch3(double(images{imnr-1}),double(images{imnr}),pt,destpos,8*threshold);
      if residual<intensitythreshold & ...
         norm(destpos-destpos2)<8*threshold,
        currpts(:,lostindex(i))=destpos2;
      end
    end
  end
  imseq{imnr}=imagedata(filenames{imnr},currpts);

  [str,mot] = bundleplc(str,mot,imseq,{'iteration=20','outputoff','lambda=1e4'});

  F12=ptobi(mot,[imnr-2,imnr-1]);F12=F12/norm(F12);
  F23=ptobi(mot,[imnr-1,imnr]);F23=F23/norm(F23);
  F13=ptobi(mot,[imnr-2,imnr]);F13=F13/norm(F13);

  %examine "unused" harris-corners
  imtmp=imagedata;
  harrispts=getpoints(imharris);
  prevpts=getpoints(imseq{imnr-1});
  for i=1:size(harrispts,2);
    [mdist,mindex]=min((harrispts(1,i)-prevpts(1,:)).^2+...
		       (harrispts(2,i)-prevpts(2,:)).^2);
    if mdist>mincornerdistance^2,
      imtmp=addpoints(imtmp,harrispts(:,i));
    end
  end


  newimseq=findmorematches(images(imnr-2:imnr),imtmp,F12,F23,F13,threshold,4*intensitythreshold/3);

  tmpmot=motion(getcameras(mot,imnr-2:imnr));
  newstr=intsecpoints(newimseq,tmpmot);
  newstr=bundleplc(newstr,tmpmot,newimseq,{'lambda=1','iteration=50','outputoff','structure'});

  %calculate reprojection error
  err=zeros(1,size(newstr,1));
  for i=1:3;
    imreproj=project(newstr,tmpmot,i);
    pts=pflat(getpoints(newimseq{i}));
    repts=pflat(getpoints(imreproj));
    err=err+sum((pts-repts).^2);
  end
  ind=find(err<6*threshold);


  for i=1:imnr;
    if i<imnr-2,
      imseq{i}=imseq{i}+imagedata([],NaN*ones(3,length(ind)));
    else
      imseq{i}=imseq{i}+imagedata([],getpoints(newimseq{i-imnr+3},ind));
    end
  end;
  str=str+structure(getpoints(newstr,ind));

  [str,mot] = bundleplc(str,mot,imseq,{'iteration=20','outputoff','lambda=1e4'});

  %FINALLY, determine which ones to keep,
  err=zeros(1,size(str,1));
  nbrappeared=zeros(1,size(str,1));
  for i=1:imnr;
    imreproj=project(str,mot,i);
    pts=pflat(getpoints(imseq{i}));
    repts=pflat(getpoints(imreproj));
    tmp=sum((pts-repts).^2);
    ind=find(~isnan(tmp));
    err(ind)=err(ind)+tmp(ind);
    nbrappeared(ind)=nbrappeared(ind)+1;
  end
  ind=find(err./nbrappeared<2*threshold);
  for i=1:imnr;
    imseq{i}=imagedata(filenames{i},getpoints(imseq{i},ind));
  end;
  str=structure(getpoints(str,ind));

  [str,mot] = bundleplc(str,mot,imseq,{'iteration=10','outputoff','lambda=1e4'});

  nbrpoints=size(imseq{1},1);


%rmspoints(str,mot,imseq);


  if length(logfile)>0,
    eval(['save ',logfile,num2str(imnr),' str mot imseq nbrpoints']);
  end;

end; %image-loop

%ii=1;figure(ii);reproject(str,mot,imseq,ii,'numbered');
