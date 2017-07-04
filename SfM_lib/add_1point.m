function [u,U] = add_1point(settings,P,U,u, imindex)

img_path = settings.img_path;
imnames = settings.imnames;

KK = settings.KK;
kc = settings.kc;
fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);

nbr = u.pointnr;

for i=imindex;
    figure(i);clf;
    filename = strcat(img_path,imnames(i).name);
    if isfield(imnames(i),'ts'),
        im = LoadImage(filename, imnames(i).ts, settings.LUT);
    else
        im = imread(filename);
    end
    imagesc(im);
    hold on;zoom on;
end
disp('zoom in');
pause;
xlist=zeros(3,0);
for i=imindex;
    figure(i);zoom off;
    x = ginput2(1);
    plot(x(1),x(2),'r*');
    x_norm = normalize(x',fc,cc,kc,alpha_c);
    xlist(:,end+1) = [x_norm;1];
    u.points{i}(:,end+1)=[x_norm;1];
    u.index{i}(:,end+1)=nbr+1;
end
u.pointnr = nbr+1;
Utmp=linf_triangulation(xlist(1:2,:),P(imindex),20,0.01);
%Utmp = intsec2views(P{imindex(1)},P{imindex(2)},xlist(:,1),xlist(:,2));
U(:,end+1)=Utmp;
for i=imindex;
    figure(i);
    utmp = P{i}*Utmp;
    utmp = utmp / utmp(3);
    utmp = KK*pextend(apply_distortion(utmp(1:2,:),kc));
    plot(utmp(1,:),utmp(2,:),'o');
end

