function [str,mot,imseq,imfiles,rgb]=visualize_import(bib,filename)
% function [str,mot,imseq,imfiles,rgb]=visualize_import(bib,filename)
% Import results from Visualize
% INPUT:
%         bib - directory
%         filename - filename of imported file (without directory).
%                    Should always have extension '.out'
%                    Hence, full filename is given by [bib,filename]
% OUTPUT:
%         str - structure
%         mot - motion
%         imfiles - filenames of images
%         rgb - 3xn matrix for point colors where n is the number of points
%


%input file
fid = fopen([bib,filename]);
nbrcameras=fscanf(fid,'%d');
imfiles=cell(1,nbrcameras);
imseq=cell(1,nbrcameras);
mot=motion;
for ii=1:nbrcameras,
    imfiles{ii}=fgetl(fid);
    imseq{ii}=imagedata([bib,imfiles{ii}]);
    
    tline=fgetl(fid);
    r1=sscanf(tline,'%f');
    tline=fgetl(fid);
    r2=sscanf(tline,'%f');
    tline=fgetl(fid);
    r3=sscanf(tline,'%f');
    tline=fgetl(fid);
    t=sscanf(tline,'%f');
    tline=fgetl(fid);
    tmp=sscanf(tline,'%f');
    focal=tmp(1);
    pp=tmp(2:3);
    K=eye(3);
    K(1,1)=focal;
    K(2,2)=focal;
    K(1:2,3)=pp;
    P=K*[[r1,r2,r3]',t];
    
    mot=addcameras(mot,P);
    
    %Closest cameras
    tline=fgetl(fid);
end
%empty line
tline=fgetl(fid);
tline=fgetl(fid);
nbrpoints=sscanf(tline,'%d');
X=zeros(3,nbrpoints);
rgb=zeros(3,nbrpoints);
for ii=1:nbrpoints,
    tline=fgetl(fid);
    X(:,ii)=sscanf(tline,'%f');
    tline=fgetl(fid);
    rgb(:,ii)=sscanf(tline,'%d');
end
str=structure(X);
fclose(fid);


