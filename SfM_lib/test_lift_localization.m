load ../data/oxford/result_shortseq1-4.mat

imnames = settings.imnames;
img_path = settings.img_path;

KK = settings.KK;
kc = settings.kc;
fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);
pixtol = 5;
%pixtol = 8

lift_threshold = 5;

test_index = [893:3187];

test_index = test_index(18)

iteration = 10000;

P0_ransac = {};

nbr_inliers_list = zeros(1,0);
dist_list = zeros(1,0);

models=20;
models=1

for kk = 1:length(test_index);
    qqq = test_index(kk);
    
    %Load test image (query)
    filename = strcat(img_path,imnames(qqq).name);
    im_orig = LoadImage(filename, imnames(qqq).ts, settings.LUT);

    load([filename num2str(imnames(qqq).ts) '.mat']);
    locs1 = kp.locs;
    descriptors1 = kp.desc;

    best_inliers = 0;
    best_P=[];best_index=[];
    
    %remove sift points in bad areas (for example, on the car)
    badindex = [];
    for ii=1:size(settings.forbidden,2),
        tmp = settings.forbidden{ii};
        badindex = union(badindex,find(locs1(1,:) >= tmp(1,1) & locs1(1,:) <= tmp(1,2) & ...
                                       locs1(2,:) >= tmp(2,1) & locs1(2,:) <= tmp(2,2)));
    end
    
    if 1,
        figure(1);clf;
        plot(imagedata(im_orig,locs1(1:2,:)));hold on;
        plot(imagedata([],locs1(1:2,badindex)),'kx');
        hold off;
        drawnow;
    end
    locs1(:,badindex)=[];
    descriptors1(:,badindex)=[];
    
    %compute normalized image coordinates
    fc = KK([1 5]);
    cc = KK(1:2,3);
    alpha_c = KK(1,2)/fc(1);
    locs1_normalized = normalize(locs1(1:2,:),fc,cc,kc,alpha_c);
    
    

    for run=1:models;
%        [qqq, run]
        eval(['load ../data/oxford/lift_models/shortseq1_model',num2str(run)]);



        %%%%%%%%%%%%%%%%%%%%
        % Create kd-tree
        %%%%%%%%%%%%%%%%%%%%
        nbr_images = length(P);
        sift = zeros(128,0);
        index3d = zeros(1,0);
        for ii=1:nbr_images,
            nbr_points = size(lift_descriptors{ii},2);
            sift(:,end+[1:nbr_points])=double(lift_descriptors{ii});
            index3d(1,end+[1:nbr_points])=u.index{ii}(1:nbr_points); %gives index of which 3D point, given a sift
        end

        %Build kdtree
        kdtree = vl_kdtreebuild(sift,'NumTrees', 12);







        %query sifts from test image
        nbr_neighbours = 1;

        qsift = descriptors1;
        [siftindex,dist] = vl_kdtreequery(kdtree,sift,qsift,'MAXCOMPARISONS',1000,'NUMNEIGHBORS',nbr_neighbours);


        %construct matches
        pp2d = zeros(2,0);
        pp3d = zeros(3,0);

        for ii=1:size(qsift,2);
            ind = unique(index3d(siftindex(:,ii)));
            for jj=1:length(ind);
                if dist(:,ii)<=lift_threshold,
                    pp2d(:,end+1)=locs1_normalized(1:2,ii);
                    pp3d(:,end+1)=U(1:3,ind(jj));
                end
            end
        end

        pp3d_ext = pextend(pp3d);

        %RANSAC for 2D-3D absolute pose
        threshold =pixtol/KK(1,1);
        nbr_matches = size(pp2d,2);

        improved=0;
        for iii=1:iteration,
         per=randperm(nbr_matches);per = per(1:3);
         pp2dselect = pp2d(:,per);
         if norm(pp2dselect(:,1)-pp2dselect(:,2))>1e-3 && ...
           norm(pp2dselect(:,1)-pp2dselect(:,3))>1e-3 && ...
           norm(pp2dselect(:,2)-pp2dselect(:,3))>1e-3,

            motlist = resec3pts(imagedata([],pp2dselect),structure(pp3d(:,per)));
            P_ransac = getcameras(motlist);
            for jj=1:length(P_ransac);
                tmp = P_ransac{jj} * pp3d_ext;
                posdepth = tmp(3,:)>0;
                tmp = pflat(tmp);
                inliers = ((sum((tmp(1:2,:)-pp2d).^2))<=threshold^2) & posdepth;
                nbr_inliers = sum(inliers);
                if nbr_inliers>best_inliers,
                    best_inliers = nbr_inliers;
                    best_P = P_ransac{jj};
                    best_index = find(inliers);
                    improved=1;
                end
            end
         end
        end
        if improved %for this model?

            %bundle the solution
            mot = motion(best_P);
            imseq = {imagedata([],pp2d(:,best_index))};
            str = structure(pp3d(:,best_index));
            [~,mot]=bundleplc(str,mot,imseq,{'motion','autocalib=11111','outputoff'});
            best_P = getcameras(mot,1);

            P0_ransac{qqq}=best_P;
            nbr_inliers_list(qqq)=best_inliers;

            C_ransac = pflat(null(best_P));
            C_gt = pflat(null(P_uncalib{qqq}));
            
            dist_list(qqq) = norm(C_gt-C_ransac);

            disp(['Image: ',num2str(qqq),', model: ',num2str(run),'. Inliers: ',num2str(best_inliers),'. Distance error of centres: ',num2str(norm(C_gt-C_ransac)),'.']);

            if 1,
                figure(3);clf;
                plot(imagedata(im_orig));hold on;
        
                tmp = pp2d(:,best_index);
                tmp = KK*pextend(apply_distortion(tmp(1:2,:),kc));
                plot(tmp(1,:),tmp(2,:),'*');
                PU = pflat(best_P*pextend(pp3d(:,best_index)));
                PU = KK*pextend(apply_distortion(PU(1:2,:),kc));
                plot(PU(1,:),PU(2,:),'ro');
                hold off;
            end

            figure(2);clf;plot(focalpoints(motion(P0_ransac(test_index(1:kk)))),'ro');hold on;plot(focalpoints(motion(P_uncalib(test_index(1:kk)))));axis equal;
        else
            disp(['Image: ',num2str(qqq),', model: ',num2str(run),'. No improvement.']);
        end

        drawnow;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    disp(['Image ',num2str(qqq),' done. Good: ',num2str(sum(dist_list(test_index(1:kk))<5e-4)),' out of ',num2str(kk)]);


end

%eval(['save ../data/oxford/lift_models/result_localization_shortseq1 dist_list nbr_inliers_list P0_ransac models iteration test_index']);
