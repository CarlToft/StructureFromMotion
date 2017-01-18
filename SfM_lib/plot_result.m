function plot_result(settings,P,U,u, imindex)

if nargin<5,
    imindex = 1:length(P);
end
img_path = settings.img_path;
imnames = settings.imnames;

KK = settings.KK;
kc = settings.kc;

figure(1);
plot3(U(1,:),U(2,:),U(3,:),'.')
axis equal;
for i=imindex;
    figure(2+mod(i-1,20));
    filename = strcat(img_path,imnames(i).name);
    if isfield(imnames(i),'ts'),
        im = LoadImage(filename, imnames(i).ts, settings.LUT);
    else
        im = imread(filename);
    end
    imagesc(im);
    hold on;
    tmp = u.points{i}; 
    tmp = KK*pextend(apply_distortion(tmp(1:2,:),kc));
    plot(tmp(1,:),tmp(2,:),'*');
    PU = pflat(P{i}*U(:,u.index{i}));
    PU = KK*pextend(apply_distortion(PU(1:2,:),kc));
    plot(PU(1,:),PU(2,:),'ro');
    ax=axis;dev=0.003*ax(2);
    for pp=1:size(PU,2);
      h=text(double(PU(1,pp)-dev),double(PU(2,pp)-dev),num2str(u.index{i}(pp)));
      set(h,'Fontsize',11);set(h,'Color','g');
    end
    hold off;
    disp('Press a key to see next image');
    pause;
end

function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);
