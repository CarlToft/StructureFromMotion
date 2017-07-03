function [curve]=add_1planarcurve(settings,u,U,P,imnr);

seg_path = settings.seg_path;

KK = settings.KK;
kc = settings.kc;
fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);


%%%%%%%%%%%%%%%%%%%%%% 
    i = imnr;

    filename = strcat(settings.img_path,settings.imnames(i).name);
    if isfield(settings.imnames(i),'ts'),
        im = LoadImage(filename, settings.imnames(i).ts, settings.LUT);
    else
        im = imread(filename);
    end

    im_seg = imread([seg_path,'im',num2str(i),'.png']);
    figure(1);clf;subplot(1,2,1);image(im);hold on;
    subplot(1,2,2);imagesc(im_seg);hold on;
    
    [m,n]=size(im_seg);
    
%    xx=pextend(ginput()');

    [~,x,y] = roipoly;
    xx=[x,y]';xx=xx(:,1:end-1);
    
    pp_im = KK*pextend(apply_distortion(u.points{i}(1:2,:),kc));
    
    [in,on] = inpolygon(pp_im(1,:),pp_im(2,:),xx(1,:),xx(2,:));
    
    plot(imagedata([],[xx,xx(:,1)]),'c-*');
    plot(imagedata([],pp_im(:,in)),'r.');
    
    thr = 0.0001;
    index = find(in);
    [pl, U_inliers, index_best] = ransac_plane(U(:,u.index{i}(index)),thr);
    index = index(index_best);
    plot(imagedata([],pp_im(:,index)),'go');
    
    if size(U_inliers,2)<8,
        error('too few');
    end
    
    ll_norm = pextend(normalize(xx,fc,cc,kc,alpha_c));
    
    lambda = (pl(1:3)'*P{i}(:,1:3)'*P{i}(:,4)-pl(4))./(pl(1:3)'*P{i}(:,1:3)'*ll_norm);
    curve = pextend(P{i}(:,1:3)'*ll_norm.*(ones(3,1)*lambda)-P{i}(:,1:3)'*P{i}(:,4)*ones(1,size(ll_norm,2)));
    
    
    

    