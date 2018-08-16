function [imseq,str] = annotate_3dpoints(settings,mot,imindex, nbrpoints, imseq, str);

nbrimages = size(settings.imnames,2);
if nargin<5,
    imseq = cell(1,nbrimages);
    for ii=1:nbrimages,
        imseq{ii}=imagedata;
    end
    str = structure;
end

img_path = settings.img_path;
imnames = settings.imnames;
KK = settings.KK;
kc = settings.kc;

%KK = [857.48296979 0 968.06224829;0 876.71824265 556.37145899;0 0 1]
%kc = [-2.57614020e-1;8.77086999e-2;-2.56970803e-4;-5.93390389e-4;-1.52194091e-2]


fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);


for ii=imindex
    figure(ii);clf;
    filename = strcat(img_path,imnames(ii).name);
    if isfield(imnames(ii),'ts'),
        im = LoadImage(filename, imnames(ii).ts, settings.LUT);
    else
        im = imread(filename);
    end
    imagesc(im);
    hold on;zoom on;
end

disp('zoom in');
pause;
for qq = 1:nbrpoints,
    xlist=NaN*ones(3,nbrimages);
    xlist2 = zeros(3,0);
    for i=imindex;
        figure(i);zoom off;
        x = ginput2(1);
        plot(x(1),x(2),'r*');
        x_norm = normalize(x',fc,cc,kc,alpha_c);
        xlist(:,i) = [x_norm;1];
        xlist2(:,end+1)=  [x_norm;1];
    end
    for i = 1:nbrimages,
        imseq{i}=addpoints(imseq{i},xlist(:,i));
    end
    Utmp=linf_triangulation(xlist2(1:2,:),getcameras(mot,imindex),20,0.01);
    str = addpoints(str,Utmp);
    for i=imindex;
        figure(i);
        utmp = getcameras(mot,i)*Utmp;
        utmp = utmp / utmp(3);
        utmp = KK*pextend(apply_distortion(utmp(1:2,:),kc));
        plot(utmp(1,:),utmp(2,:),'o');
    end
end
    
    
    
end
