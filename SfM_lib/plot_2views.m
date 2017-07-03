function plot_2views(settings,P,U,u, imindex, click)

if nargin<6,
    click=0;
end
if nargin<5,
    imindex = 1:2;
end
img_path = settings.img_path;
imnames = settings.imnames;

KK = settings.KK;
kc = settings.kc;

pts_index = intersect(u.index{imindex(1)},u.index{imindex(2)});

for i=imindex;
    figure(i);clf;
    filename = strcat(img_path,imnames(i).name);
    if isfield(imnames(i),'ts'),
        im = LoadImage(filename, imnames(i).ts, settings.LUT);
    else
        im = imread(filename);
    end
    imagesc(im);
    hold on;
    [the_index,ind]=intersect(u.index{i},pts_index);
    tmp = u.points{i}(:,ind); 
    tmp = KK*pextend(apply_distortion(tmp(1:2,:),kc));
    plot(tmp(1,:),tmp(2,:),'*');
    unflat = P{i}*U(:,the_index);
    depths = unflat(3,:);
    PU = pflat(unflat);
    PU = KK*pextend(apply_distortion(PU(1:2,:),kc));
    plot(PU(1,:),PU(2,:),'ro');
    ax=axis;dev=0.003*ax(2);
    for pp=1:size(PU,2);
      h=text(double(PU(1,pp)-dev),double(PU(2,pp)-dev),num2str(the_index(pp)));
      set(h,'Fontsize',11);set(h,'Color','g');
    end
    axis equal;
end

for ii=1:click,
    fc = KK([1 5]);
    cc = KK(1:2,3);
    alpha_c = KK(1,2)/fc(1);
    figure(imindex(1));
    x = ginput2(1);
    plot(x(1),x(2),'r*');
    
    x_norm = normalize(x',fc,cc,kc,alpha_c);
    F=ptobi(motion(P(imindex(1:2))));
    F=F/norm(F);
    ll= F'*pextend(x_norm);
    
    x_other=zeros(3,0);
    for xx=-2:0.01:2
        x_other(:,end+1) = pflat(cross(ll,[1,0,xx]'));
    end
    
    PU = KK*pextend(apply_distortion(x_other(1:2,:),kc));
    figure(imindex(2));ax=axis;
    plot(PU(1,:),PU(2,:),'r-');
    axis(ax);
end


function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);
