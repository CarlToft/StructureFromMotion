function annotate_plot(settings,P_uncalib,imindex, tri, Utri, triobj, Utriobj,labelobj);

nbrimages = size(settings.imnames,2);

img_path = settings.img_path;
imnames = settings.imnames;
KK = settings.KK;
kc = settings.kc;
fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);

cols = 'rgbyckw';

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
    
    if size(tri,2)>0,
        pp = pflat(P_uncalib{ii}*pextend(Utri));
        pp = KK*pextend(apply_distortion(pp(1:2,:),kc));
        plot(imagedata([],pp),'bx');
        for jj = 1:size(tri,2);
            index = [tri(:,jj);tri(1,jj)];
            plot(pp(1,index),pp(2,index),'y-');
        end
    end
    if size(triobj,2)>0,
        pp = pflat(P_uncalib{ii}*pextend(Utriobj));
        pp = KK*pextend(apply_distortion(pp(1:2,:),kc));
        plot(imagedata([],pp),'r.');
        for jj = size(triobj,2):-1:1;
            index = [triobj(:,jj);triobj(1,jj)];
            plot(pp(1,index),pp(2,index),[cols(mod(labelobj(jj),length(cols))+1),'-']);
        end
    end
end
