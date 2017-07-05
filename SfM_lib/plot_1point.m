function plot_1point(settings,P,U,u, Uindex, always_plot_in_image)

if nargin<6,
    always_plot_in_image = 0;
end
    
img_path = settings.img_path;
imnames = settings.imnames;

KK = settings.KK;
kc = settings.kc;


for i=1:length(P);
    ind=find(u.index{i}==Uindex);
    if ~isempty(ind) || always_plot_in_image,
        figure(i);clf;
        filename = strcat(img_path,imnames(i).name);
        if isfield(imnames(i),'ts'),
            im = LoadImage(filename, imnames(i).ts, settings.LUT);
        else
            im = imread(filename);
        end
        imagesc(im);
        hold on;
        if ~isempty(ind)
            tmp = u.points{i}(:,ind); 
            tmp = KK*pextend(apply_distortion(tmp(1:2,:),kc));
            plot(tmp(1,:),tmp(2,:),'*');
        end
        unflat = P{i}*U(:,Uindex);
        depths = unflat(3,:);
        PU = pflat(unflat);
        PU = KK*pextend(apply_distortion(PU(1:2,:),kc));
        plot(PU(1,:),PU(2,:),'ro');
        ax=axis;dev=0.003*ax(2);
        for pp=1:size(PU,2);
            h=text(double(PU(1,pp)-dev),double(PU(2,pp)-dev),num2str(Uindex(pp)));
            set(h,'Fontsize',11);set(h,'Color','g');
        end
        axis equal;
    end
end


function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);
