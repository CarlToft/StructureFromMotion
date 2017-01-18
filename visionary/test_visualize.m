
%Test program for exporting files to visualizer
%There are two cases, first second example is synthetic example,
%  the first example shows how to import from bundler and then export

% CASE 1: import a scene from bundler and then export

bibdir='C:\bundler\bundler-v0.3-binary\examples\ET';
[str,mot,Kmot,f,pp,kk, imseq0, imundist, rgb]=bundler_import(bibdir);

nbrcameras=size(mot);
imfiles=cell(1,nbrcameras);

for ii=1:nbrcameras,
    tmp=getfilename(imseq0{ii});
    imfiles{ii}=tmp(length(bibdir)+2:end);
end
visualize_export([bibdir,'\ET.out'],str,Kmot,imfiles,rgb);

% CASE 2: generate a synthetic scene as well as images and then export
filename='synt.out'; %export filename
imagename='syntim';  %image filename

nbrcameras=10;
nbrpoints=50; %random points
background=0.8; %background color
sx=300; %size of images
sy=200;

rgb=rand(3,nbrpoints);
pp=[sx/2;sy/2]; %principal point
margin=2; %size of projected point in image

[str,mot,imseq]=randomscene(nbrcameras,nbrpoints,0,0,{'noise=0','zeroprincipal'});
T=diag([1/100,1/100,1/100,1]);
mot=changecsystem(mot,T);
str=changecsystem(str,T);

T=diag([1/5,1/5,1]);
T(1:2,3)=pp;
mot=changecsystem(mot,T);
imnames=cell(1,nbrcameras);
for ii=1:nbrcameras,
    imseq{ii}=changecsystem(imseq{ii},T);
    im=background*ones(sy,sx,3);
    
    pt=round(pflat(project(str,mot,ii)));
    index=find(pt(1,:)>margin & pt(1,:)<sx-margin & pt(2,:)>margin & pt(2,:)<sy-margin);
    for jj=index,
        im(pt(2,jj)+[-margin:margin],pt(1,jj)+[-margin:margin],1)=rgb(1,jj);
        im(pt(2,jj)+[-margin:margin],pt(1,jj)+[-margin:margin],2)=rgb(2,jj);
        im(pt(2,jj)+[-margin:margin],pt(1,jj)+[-margin:margin],3)=rgb(3,jj);
    end
    imseq{ii}=imagedata(im,getpoints(imseq{ii}));
    
    imfile=sprintf('./%s%03d.jpg',imagename,ii);
    imwrite(im,imfile);
    imfiles{ii}=imfile;
end

visualize_export(filename,str,mot,imfiles,rgb);
