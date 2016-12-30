function [i,crf]=harris(im,option);
% [i,crf]=HARRIS(im,option) - creates an IMAGEDATA object with corners using HARRIS
% im: image matrix
% option: (cell of strings)
% - 'threshold=X', sets corner response threshold to X
% - 'nbrcorners=X', gets the X strongest corners in the image
% - 'noimage', image is not included in IMAGEDATA object i
% - 'mindistance=X' - minimum pixel distance between corners
% Output:
%  i - imagedata object with corner points
%  crf - corner response function
% Default is 'threshold=1', 'mindistance=20'
% See also: harrispro.m

cornerthreshold=1;
nbrcorners=0;
includeimage=1;

mindistance=20;

if nargin>=2,
  if strmatch('threshold=',option);
    q=strmatch('threshold=',option);
    strthreshold = option{q}(11:length(option{q}));
    cornerthreshold=str2num(strthreshold);
  end
  if strmatch('nbrcorners=',option);
    q=strmatch('nbrcorners=',option);
    strnbrcorners = option{q}(12:length(option{q}));
    nbrcorners=str2num(strnbrcorners);
  end
  if strmatch('mindistance=',option);
    q=strmatch('mindistance=',option);
    mindistance=str2num(option{q}(13:length(option{q})));
  end
  if strmatch('noimage',option);
    includeimage=0;
  end
end %option

a1=1.5;
filtsize1=round(2*a1);
[x,y]=meshgrid(-filtsize1:filtsize1,-filtsize1:filtsize1);
filtx=-2*x.*exp(- (x.^2 + y.^2)/a1^2)/(a1^4*pi);
filty=-2*y.*exp(- (x.^2 + y.^2)/a1^2)/(a1^4*pi);
Lx=conv2(im,filtx,'valid');
Ly=conv2(im,filty,'valid');

a2=1.5;
filtsize2=round(2*a2);
[x,y]=meshgrid(-filtsize2:filtsize2,-filtsize2:filtsize2);
filt=exp(- (x.^2 + y.^2)/a2^2)/(a2^2*pi);

Wxx=conv2(Lx.*Lx,filt,'valid');
Wxy=conv2(Lx.*Ly,filt,'valid');
Wyy=conv2(Ly.*Ly,filt,'valid');

%Wdet=Wxx.*Wyy-Wxy.^2;
%Wtr=Wxx+Wyy;

ignoreborder=10;
%utbild=abs(Wxx.*Wyy-Wxy.^2);
utbild=abs((Wxx.*Wyy-Wxy.^2)./(Wxx+Wyy+eps));
utbild=utbild(1+ignoreborder:end-ignoreborder,1+ignoreborder:end-ignoreborder);


% Lite extra test
%[m,n]=size(utbild);
%mask=zeros(m,n);
%cutsize=((filtsize1+filtsize2)*2+1);
%mask(cutsize:m-cutsize,cutsize:n-cutsize)=ones(m-2*cutsize+1,n-2*cutsize+1);
%utbild = mask.*utbild;

cont=1;
topplista=[];
crf=[];
[m,n]=size(utbild);

shifted=ignoreborder+filtsize1+filtsize2;
while cont,
  [xmax,ypos]=max(utbild);
  [ymax,xpos]=max(xmax);
  ymax=ymax/9; %normalize
  ypos=ypos(xpos);

  if nbrcorners==0,
    if ymax>cornerthreshold,
  	topplista=[topplista,(shifted+[xpos;ypos])];
        crf=[crf,ymax];
    else
	cont=0;
    end
  else %pick a fixed number of corners
    topplista=[topplista,(shifted+[xpos;ypos])];
    crf=[crf,ymax];
    nbrcorners=nbrcorners-1;
    if nbrcorners==0,
	cont=0;
    end
  end
  ind1=max(ypos-mindistance,1):min(ypos+mindistance,m);
  ind2=max(xpos-mindistance,1):min(xpos+mindistance,n);

  utbild(ind1,ind2)=zeros(length(ind1),length(ind2));
end;

if includeimage,
  i=imagedata(im,topplista);
else
  i=imagedata([],topplista);
end
