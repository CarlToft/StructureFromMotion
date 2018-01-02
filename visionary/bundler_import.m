function [str, mot, Kmot, f, pp, kk, im, imundist, rgb] = bundler_import(bibdir)
%Import results from Noah Snavely's Bundler (Photo Tourism) to visionary
% INPUT: filename
% OUTPUT: str - structure
%         mot - motion (without Calibration matrix)
%         Kmot - motion (with calibration matrix)
%         f - focal lengths
%         pp - principal points
%         kk - radial distortions
%         im - cell of imagadata, original coords & image files
%         imundist - cell of imagedata, undistorted, normalized coords
%         rgb - 3xn matrix of point colors
%
% See "test_bundler" for examples
%
filename = [bibdir,'/bundle/bundle.out'];
fid = fopen(filename);

fscanf(fid, '%s', 4);

nbr_cameras = fscanf(fid, '%d',1);
nbr_points = fscanf(fid, '%d',1);

% Set up camera output
mot = motion;

f=zeros(1,0);
kk=zeros(2,0);
index=zeros(1,nbr_cameras);
im={};
im0={};
for camera = 1:nbr_cameras
    focal_length = fscanf(fid, '%f',1);
    k1k2 = fscanf(fid, '%f', 2);
    R = fscanf(fid, '%f', [3,3])';
    t = fscanf(fid, '%f', 3);
    if focal_length ~=0
%        mot = mot + motion(-[focal_length 0 0; 0 focal_length 0; 0 0 1]*[R t]);
        mot = mot + motion(diag([1,-1,-1])*[R t]);
        f(end+1)=focal_length;
        kk(:,end+1)=k1k2';
        index(camera)=length(f);
        im0{end+1}=NaN*ones(2,nbr_points);
    end
end

nbr_actual_cameras=size(mot);

str = structure;
rgb=zeros(3,nbr_points);

for point_ind = 1:nbr_points
    % Read the position
    point = fscanf(fid, '%f', 3);
    
    str = str + structure(point);
    
    % rgb color
    rgbtmp = fscanf(fid, '%f', 3);
    rgb(:,point_ind)=rgbtmp;
    
    % 2d information
    nbr_views = fscanf(fid, '%d', 1);
    tmp=reshape(fscanf(fid, '%f', 4*nbr_views),4,nbr_views);
    for ii=1:nbr_views,
        imindex=index(tmp(1,ii)+1);
        im0{imindex}(:,point_ind)=diag([1,-1])*tmp(3:4,ii);
    end
end

% Close file
fclose(fid);



% Open "list.txt"
filename = [bibdir,filesep,'list.txt'];
fid = fopen(filename);

imnames={};
for ii=1:nbr_cameras,
    % Read image filename
    
    tline = fgetl(fid);
    imnames{ii}=sscanf(tline, '%s', 1);
end
% Close file
fclose(fid);


for ii=1:length(im0),
    tmpindex=find(index==ii);
    tmpfile = [bibdir,filesep,imnames{tmpindex}];
    tmpim=imread(tmpfile);
    pp(1,ii)=size(tmpim,2)/2;
    pp(2,ii)=size(tmpim,1)/2;
    tmp=im0{ii};
    tmp(1,:)=tmp(1,:)+pp(1,ii);
    tmp(2,:)=tmp(2,:)+pp(2,ii)-1;
    im{ii}=imagedata(tmpfile,tmp);
end

imundist=undist_radial(im0,f,kk);

Kmot=motion;
for ii=1:size(mot),
    Puncalib=getcameras(mot,ii);
    K=diag([f(ii),f(ii),1]);
    K(1:2,3)=pp(:,ii);
    Kmot=addcameras(Kmot,K*Puncalib);
end

%normalize coordinate system
X=pflat(str);X(4,:)=[];tr=mean(X')';
sc=1/std(X(:));
T=sc*eye(4);T(4,4)=1;
T(1:3,4)=-sc*tr;
str=changecsystem(str,T);
mot=changecsystem(mot,T);
Kmot=changecsystem(Kmot,T);


%imdist=dist_radial(imundist,f,kk)
%ii=13;imindex=2;tmp=pflat(P{imindex}*X(:,ii));tmp(3)=[];l=norm(tmp);[tmp*(1+l^2*kk(1,imindex)+l^4*kk(2,imindex))*f(imindex),im0{imindex}(:,ii)]
%ii=13;imindex=2;tmp=pflat(P{imindex}*X(:,ii));[tmp,getpoints(imundist{imindex},ii)]
%1+l^2*kk(1,imindex)+l^4*kk(2,imindex))*f(imindex),im0{imindex}(:,ii)]
