function [s,m,residual,resnew,lambda] = bundleplc(s,m,imseq,option);
% function [s,m,redidual,resnew] = bundleplc(s,m,imseq,option);
% Bundle adjustment for combinations of points, lines and conics
% Input:
%   s - structure object
%   m - motion object
%   imseq - cell list of imagedata objects
%   option: (cell of strings)
%     'affine' - affine s & m
%     'iteration=X' - iterate X times
%     'lambda=X' - set lambda to X
%     'outputoff' - No results displayed
%     'autocalib=XXXXX' - Autocalibration
%       specify for each of focal length, aspect ratio, skew,
%       principal point x, principal point y,
%       if the parameter is
%        1 known and nominal,
%        2 unknown but constant in the sequence,
%        3 unknown and varying.
%     'structure' - Only optimize structure
%     'motion' - Only optimize motion
%     'knownrotation' - Only optimize camera centres and structure
%     'focal=[f,sigma]' - Focal length at f with std sigma
%     'ar=[a,sigma]' - Aspect ratio at a with std sigma
%     'pp=[u,v,sigma]' - Principal point at (u,v) with std sigma
%     'lockcsystem' - Lock coordinate system. Locks first camera
%        and one of second camera centre coordinates. NB: Works
%        only in autocalibration mode.
%     'pflat' - 3D points are pflat'ed (default is psphere).
%     'stdlines=X' - standard deviation of line noise
%     'zplane=[index1,index2,...] - index of points on z=0 plane
% Output:
%   s        - updated structure
%   m        - updated motion
%   residual - total reprojection error
%   resnew   - calculated residuals
%   Default is projective bundle.
warningstate=warning('off');

if ~isa(s,'structure'),
    s=structure(s);
end
if ~isa(m,'motion'),
    m=motion(m);
end

caliboptions=[];
mode = 1;
modestr={'projective','affine','autocalibration','structure','motion','autocalibration+motion','knownrotation'};
iteration=10;
startlambda=10;
outputon=1;
focal=[];
ar=[];
pp=[];
lockcsystem=0;
lockstr={'coordinate system free','coordinate system locked'};
pflatpoints=0;
stdlines=1; % standard deviation for lines (endpoint distance)

zplane=[];
if nargin>3,
  if strmatch('affine',option),
    mode = 2;
  end
  if strmatch('outputoff',option),
    outputon=0;
  end
  if strmatch('iteration=',option),
    q=strmatch('iteration=',option);
    striter = option{q}(11:length(option{q}));
    iteration=str2num(striter);
  end
  if strmatch('lambda=',option),
    q=strmatch('lambda=',option);
    strlambda = option{q}(8:length(option{q}));
    startlambda=str2num(strlambda);
  end
  if strmatch('autocalib=',option),
    mode = 3;
    q=strmatch('autocalib=',option);
    strautocalib = option{q}(11:length(option{q}));
    caliboptions=[str2num(strautocalib(1)) str2num(strautocalib(2)) str2num(strautocalib(3)) str2num(strautocalib(4)) str2num(strautocalib(5))];
  end
  if strmatch('structure',option),
    mode = 4;
  end
  if strmatch('motion',option),
    if mode==3,
      mode=6;   %autocalibration+motion
    else
      mode = 5; %motion (i.e. projective)
    end
  end
  if strmatch('knownrotation',option),
      mode=7;
  end
  if strmatch('focal=',option) & ~isempty(caliboptions),
    if caliboptions(1)>1,
      q=strmatch('focal=',option);
      tmp=option{q};
      tmp(end+1)=';';
      eval(tmp);
      focal(3)=caliboptions(1);
    end
  end
  if strmatch('ar=',option) & ~isempty(caliboptions),
    if caliboptions(2)>1,
      q=strmatch('ar=',option);
      tmp=option{q};
      tmp(end+1)=';';
      eval(tmp);
      ar(3)=caliboptions(2);
    end
  end
  if strmatch('pp=',option) & ~isempty(caliboptions),
    if sum(caliboptions(4:5)>1)>0,
      q=strmatch('pp=',option);
      tmp=option{q};
      tmp(end+1)=';';
      eval(tmp);
      pp(4:5)=caliboptions(4:5);
    end
  end
  if strmatch('lockcsystem',option) & mode == 3,
    lockcsystem=1;
  end
  if strmatch('pflat',option),
    pflatpoints=1;
  end
  if strmatch('stdlines=',option),
    q=strmatch('stdlines=',option);
    'I accidently removed these lines....'
    'redo them before use....'
    keyboard;
  end
  if strmatch('zplane=',option),
    q=strmatch('zplane=',option);
    tmp=option{q};
    tmp(end+1)=';';
    eval(tmp);
  end
end
% Extract points, lines and conics from imseq
antalbilder = length(imseq);
fs=size(imseq{1});
nbrpts=fs(1); nbrlines=fs(2); nbrconics=fs(3);
pointindex=zeros(1,nbrpts);
lineindex=zeros(1,nbrlines);
conicindex=zeros(1,nbrconics);

nbrfeatures = nbrpts + nbrlines + nbrconics;
mm=cell(antalbilder,1);
ss=cell(1,0);
im=cell(antalbilder,0);

if mode==3 | mode==6,
%  P1=getcameras(m,1);
%  K1=rq(P1);
%  T=[inv(K1)*P1;0,0,0,1];
%  m=changecsystem(m,T);
%  s=changecsystem(s,T);

  for i=1:antalbilder,
    Pi=getcameras(m,i);
    Ki=rq(Pi);
    mm{i}=pmatrix(Pi/Ki(3,3));
  end
else
  for i=1:antalbilder,
    temp=getcameras(m,i);if mode~=2, temp=temp/norm(temp);end
    mm{i}=pmatrix(temp);
  end
end

nbrpts0=0;nbrlines0=0;nbrconics0=0;
% Add points
for j=1:nbrpts,
 temp=getpoints(s,j);
 if ~isempty(find(j==zplane))>0,
    temp(3)=0;
 end
 %check that it is visible
 visible=0;
 for i=1:antalbilder,
   sl=getpoints(imseq{i},j);
   visible=max(visible,~isnan(sl(1)));
 end
 if ~isnan(temp(1)) & visible,
  nbrpts0=nbrpts0+1;
  pointindex(j)=nbrpts0;
  if pflatpoints==1,
    ss{nbrpts0}=obpoint(pflat(temp));
  else
    ss{nbrpts0}=obpoint(psphere(temp));
  end
  for i=1:antalbilder,
    [temp,L] = getpoints(imseq{i},j);
    if ~isnan(temp(1)),
      im{i,nbrpts0} = impoint(pflat(temp),L,[0;0;1]);
    end;
  end;
 end
end;
joffset=nbrpts0;

% Add lines
for j=1:nbrlines,
 temp=getlines(s,j);
 if ~isnan(temp(1)),
  nbrlines0=nbrlines0+1;
  lineindex(j)=nbrlines0;
  temp=pflat(reshape(temp,4,2));
  ss{nbrlines0+joffset}=obline(temp);
  for i=1:antalbilder,
%    [temp,L] = gethomogeneouslines(imseq{i},j);
    temp = gethomogeneouslines(imseq{i},j);
    if ~isnan(temp(1)),
      temp=temp/norm(temp(1:2));
%      im{i,nbrlines0+joffset} = imline(temp,L,[temp(1:2,1);0]);
      im{i,nbrlines0+joffset} = imline(temp,stdlines);
    end;
  end;
 end;
end;
joffset=joffset+nbrlines0;

% Add conics
for j=1:nbrconics,
 temp=getquadrics(s,j);
 if ~isnan(temp(1)),
  nbrconics0=nbrconics0+1;
  conicindex(j)=nbrconics0;
  temp=v2m(temp);
  ss{nbrconics0+joffset}=obconic(temp);
  for i=1:antalbilder,
    [temp,L] = getconics(imseq{i},j);
    if ~isnan(temp(1)),
        temp=pflat(temp);
        ntemp=[0,0,0,0,0,1]';
%      temp=temp/norm(temp);
      im{i,nbrconics0+joffset} = imconic(temp,L,ntemp);
    end;
  end;
 end
end;
cont=1;
lambda=startlambda;

if outputon==0,
  loops=iteration;
else
  loops=0;
end;

if lockcsystem>0,
  c1=pflat(null(pdp(mm{1})));c2=pflat(null(pdp(mm{2})));
  [slask,lockcsystem]=max(abs(c1(1:3)-c2(1:3)));
end


%calculate residuals
[mi,si,imi,offsetmotion,offsetsm,offsetres]=calcindexplc(im,ss,mm,mode,zplane);
resnew=calcresplc(im,ss,mm,imi,si,mi,offsetres,offsetsm,pp,focal,ar);
if outputon,
  disp(['Current error: ',num2str(norm(resnew)/sqrt(length(resnew)))]);
  plot(resnew,'.');title('Residuals');
  disp(['Mode: ',modestr{mode},', ',lockstr{(lockcsystem>0)+1}]);
end;

while cont,
%%%%%%%%%%%%%
% BUNDLE LOOOP
%%%%%%%%%%%%%
 for i=loops:-1:1;
   if outputon,
     disp(['Iterations left: ' num2str(i)]);
   end
   [ss,mm,step,lambda,oldf,newf,resnew]=bundlestepplc(im,ss,mm,lambda,mode,caliboptions,lockcsystem,pp,focal,ar,zplane);
   if outputon,
      disp(['Old error: ',num2str(oldf/sqrt(length(resnew))),'   New error: ',num2str(newf/sqrt(length(resnew)))]);
      disp(['Lognorm of gradient step: ' num2str(log(norm(step)))]);
      plot(resnew,'.');title('Residuals');drawnow;
   end
 end
%%%%%%%%%%%%%
% BUNDLE LOOOP - end
%%%%%%%%%%%%%

 if outputon,
   disp(['Current lambda: ',num2str(lambda),'   Nbr of iterations: ',num2str(iteration)]);
   askuser=1;
 else
   askuser=0;
   cont=0;
 end;

%%%%%%%%%%%%%
% ASK USER
%%%%%%%%%%%%%
 while askuser,
   disp('(c)ontinue, (s)et lambda, (i)terations, (k)eyboard');
   opt=input('(h)istogram, (l)ock/unlock coordinate system, (q)uit? ','s');
   switch opt
     case 'c', askuser=0;
     case 's', lambda=input('Lambda: ');
     case 'i', iteration=input('Iterations: ');
     case 'k', keyboard;
     case 'h', hist(resnew,21);title('Histogram of residuals');estimestd=sqrt(resnew'*resnew/(size(resnew,1)-3*nbrpts0-4*nbrlines0-9*nbrconics0-15));disp(['Estimated standard deviation: ',num2str(estimestd)]);
     case 'l', if lockcsystem==0, c1=pflat(null(pdp(mm{1})));c2=pflat(null(pdp(mm{2}))); [slask,lockcsystem]=max(abs(c1(1:3)-c2(1:3))); else lockcsystem=0;end;if outputon,disp(['Mode: ',modestr{mode},', ',lockstr{(lockcsystem>0)+1}]);end;
     case 'q', askuser=0;cont=0;
   end
 end;
% END of askuser

 loops=iteration;
end; %END of cont


% Update structure s and motion m

s=structure;
m=motion;

for j=1:nbrpts;
  if pointindex(j)~=0,
    s=addpoints(s,udu(ss{pointindex(j)}));
  else
    s=addpoints(s,NaN*ones(4,1));
  end
end;
joffset=nbrpts0;
for j=1:nbrlines;
  if lineindex(j)~=0,
    s=addlines(s,reshape(udu(ss{lineindex(j)+joffset}),8,1));
  else
    s=addlines(s,NaN*ones(8,1));
  end
end;
joffset=nbrpts0+nbrlines0;
for j=1:nbrconics;
  if conicindex(j)~=0,
    s=addquadrics(s,m2v(udu(ss{conicindex(j)+joffset})));
  else
    s=addquadrics(s,NaN*ones(10,1));
  end
end;

for i=1:antalbilder;
 m=addcameras(m,pdp(mm{i}));
end;



residual=norm(resnew);
warning(warningstate);
