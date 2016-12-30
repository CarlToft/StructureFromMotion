function visualize_export(filename,str,mot,imfiles,rgb,surface)
% function visualize_export(filename,str,mot,imfiles,rgb,surface)
% Export results to Visualize
% INPUT:
%         filename - filename of exported file.
%                    Should always have extension '.out'
%         str - structure
%         mot - motion
%         imfiles - filenames of images OR cell of imagedata objects
%         rgb - 3xn matrix for point colors where n is the number of points
%                   Should be between 0-255.
%         surface - triangulated 3D-surface
%                 - surface.vertices (nx3) containes a set of 3D-points
%                 - surface.faces (mx3) containes three indices (into the vertic list)
%                   for each triangle.
%                 - surface.facevertexcdata (mx3) containds the color (rgb)
%                   for each triangle
%
% See "test_visualize" for examples
% NB: The scale should be set so the maximum distance in 3D is around 10 units

nbrcameras=size(mot);
nbrpoints=size(str,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute closest cameras
%ORDER: left,down,forward,right,up,backward
distthreshold=1; %maximum distance for a camera to be called "close"

Pmot=getcameras(mot);
C=pflat(focalpoints(mot));C(4,:)=[];
camindex=-ones(nbrcameras,6);
for ii=1:nbrcameras,
    [K,P]=rq(Pmot{ii});
    R=P(:,1:3);
    vC=C-C(:,ii)*ones(1,nbrcameras);
    dist=sum(vC.^2);
    %determine which axis (with sign) it is closest to:
    tmp=R*vC;
    [maxc,maxindex]=max(abs(tmp));
    negindex=find(maxc>max(tmp));
    maxindex(negindex)=maxindex(negindex)+3;
    maxindex(ii)=NaN;
    
    %find closest camera in each six categories
    for jj=1:6,
        index=find(maxindex==jj);
        [slask,minindex]=min(dist(index));
        theindex=index(minindex);
        if length(theindex)==1 & sqrt(slask)<distthreshold,
            camindex(ii,jj)=theindex;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%write to file
fid = fopen(filename,'W');
fprintf(fid,'%d\n',nbrcameras);
for ii=1:nbrcameras,
    if isa(imfiles{ii},'imagedata'),
        fname=getfilename(imfiles{ii});
    else
        fname=imfiles{ii};
    end
    fprintf(fid,'%s\n',fname);
    [K,P]=rq(getcameras(mot,ii));
    K=K/K(3,3);
    focal=K(1,1);
    pptmp=K(1:2,3);
    fprintf(fid,'%f %f %f\n',P(1,1:3));
    fprintf(fid,'%f %f %f\n',P(2,1:3));
    fprintf(fid,'%f %f %f\n',P(3,1:3));
    fprintf(fid,'%f %f %f\n',P(1:3,4));
    
    %focal,px,py
    fprintf(fid,'%f %f %f\n',[focal,pptmp(:)']);
    %Closest cameras
    fprintf(fid,'%f %f %f %f %f %f\n',camindex(ii,:));
end
%empty line
fprintf(fid,'\n');
fprintf(fid,'%d\n',nbrpoints);
X=pflat(str);X(4,:)=[];
for ii=1:nbrpoints,
    fprintf(fid,'%f %f %f\n',X(:,ii));
    fprintf(fid,'%d %d %d\n',rgb(:,ii));
end

if nargin > 5
    %Triangulerad yta
    fprintf(fid,'\n');
    fprintf(fid,'%d\n', size(surface.vertices,1));
    %punkter
    for i = 1:size(surface.vertices,1);
        fprintf(fid,'%f %f %f\n',surface.vertices(i,:));
    end
    %triangelhörn index och färg
    fprintf(fid,'\n');
    fprintf(fid,'%d\n', size(surface.faces,1));
    for i = 1:size(surface.faces,1);
        fprintf(fid,'%d %d %d ',surface.faces(i,:));
        fprintf(fid,'%d %d %d\n',surface.facevertexcdata(i,:));
    end
end

fclose(fid);
