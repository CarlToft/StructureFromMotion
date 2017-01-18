function pairwise_matching(settings, seq)
imnames = settings.imnames;
img_path = settings.img_path;
scale = settings.scale;
SIFT_path = settings.SIFT_path;
save_path = settings.save_path;
KK = settings.KK;
kc = settings.kc;
distRatio = settings.distRatio; %Same rejection criterion as Lowe = 0.5

if nargin<2,
    seq = 1:length(imnames);
end
imnames = imnames(seq);
camera_graph = settings.camera_graph(seq,seq);

if ~isempty(settings.expectedEpipole),
    epipole = KK*pextend(apply_distortion(settings.expectedEpipole,kc));
end

%Compute SIFT for all the images
SIFT = cell(1,length(seq));

for i=1:length(seq);
    tmp = [i length(seq)];
    
    filename = strcat(img_path,imnames(i).name);
    if isfield(imnames(i),'ts'),
        im = LoadImage(filename, imnames(i).ts, settings.LUT);
    else
        im = imread(filename);
    end
    im = imresize(im,scale);
    
    if 0,
        [~,descriptors1,locs1] = mod_sift(im(1:1000,:,:),settings);
        [~,descriptors2,locs2] = mod_sift(im(901:end,:,:),settings);
        locs2(:,1)=locs2(:,1)+900;
    
        index = find(locs1(:,1)>900);
        for ii=index',
            index2 = find(abs(locs2(:,1)-locs1(ii,1))<1 & abs(locs2(:,2)-locs1(ii,2))<1);
            locs2(index2,:)=[];
            descriptors2(index2,:)=[];
        end
        SIFT{i}.desc = [descriptors1;descriptors2];
        SIFT{i}.locs = [locs1;locs2];
    else
        [locs1,descriptors1] = vl_sift(single(rgb2gray(im)),'PeakThresh',settings.PeakThresh,'EdgeThresh',settings.EdgeThresh);
        
        if settings.rotationvariantSIFT==1,
            %force orientations fixed
            [locs1,descriptors1] = vl_sift(single(rgb2gray(im)),'frames',[locs1(1:3,:);zeros(1,size(locs1,2))]);
        end
        SIFT{i}.desc = descriptors1';
        SIFT{i}.locs = locs1([2,1],:)';
    end
    
    SIFT{i}.locs(:,1:2) = SIFT{i}.locs(:,1:2)/scale;
    
    badindex = [];
    for ii=1:size(settings.forbidden,2),
        tmp = settings.forbidden{ii};
        badindex = union(badindex,find(SIFT{i}.locs(:,2) >= tmp(1,1) & SIFT{i}.locs(:,2) <= tmp(1,2) & ...
                                        SIFT{i}.locs(:,1) >= tmp(2,1) & SIFT{i}.locs(:,1) <= tmp(2,2)));
    end
    if ~isempty(settings.expectedEpipole),
        badindex = union(badindex,find((SIFT{i}.locs(:,2)-epipole(1)).^2+(SIFT{i}.locs(:,1)-epipole(2)).^2<=settings.epipoledistance^2));
    end
    
    if 0,
        figure(1);clf;
        p1 = SIFT{i}.locs(:,[2,1])';
        plot(imagedata(strcat(img_path,imnames(i).name),p1));hold on;
        plot(imagedata([],p1(:,badindex)),'ro');
        hold off;
    end
    SIFT{i}.locs(badindex,:)=[];
    SIFT{i}.desc(badindex,:)=[];
    
end

KK=settings.KK;
fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);
kc = settings.kc;

%Match all pairs.
pairwiseEst = cell(length(seq),length(seq));
for i = 1:length(seq);
    %Build kdtree
    tmp1 = normr(double(SIFT{i}.desc))';
    if size(tmp1,2)==0,
        for j = i+1:length(seq);
            pairwiseEst{i,j}.ind1 = [];
            pairwiseEst{i,j}.ind2 = [];
        end
    else
     kdtree = vl_kdtreebuild(tmp1,'NumTrees', 12);
     for j = i+1:length(seq);
        if ~(isempty(camera_graph)); %If there is a camera graph only use those pairs.
            if camera_graph(i,j) == 1
                [i j length(seq)]
                
                tmp2 = normr(double(SIFT{j}.desc))';
                [index,dist] = vl_kdtreequery(kdtree,tmp1,tmp2,'MAXCOMPARISONS',50,'NUMNEIGHBORS',2);
                maxvals = sum(tmp1(:,index(1,:)).*tmp2,1)';
                secondmaxvals = sum(tmp1(:,index(2,:)).*tmp2,1)';
%                maxvals = sum(SIFT{i}.desc(index(1,:),:).*SIFT{j}.desc,2);
%                secondmaxvals = sum(SIFT{i}.desc(index(2,:),:).*SIFT{j}.desc,2);
                
                matches = zeros(1,size(SIFT{j}.desc,1));
                
                tmp = acos(maxvals) < distRatio*acos(secondmaxvals);
                
%                if ~isempty(settings.expectedF);
%                    F = settings.expectedF;
%                    p1 = pextend(normalize(SIFT{i}.locs(index(1,:),[2,1])',fc,cc,kc,alpha_c));
%                    p2 = pextend(normalize(SIFT{j}.locs(:,[2,1])',fc,cc,kc,alpha_c));
%                    l1 = p1'*F;l1norm=sqrt(l1(:,1).^2+l1(:,2).^2);d1 = abs(diag(l1 * p2)./l1norm);
%                    l2 = p2'*F';l2norm=sqrt(l2(:,1).^2+l2(:,2).^2);d2 = abs(diag(l2 * p1)./l2norm);
%                    dist=d1+d2;
%                    tmp = tmp & (dist < 5*settings.RANSAC_pixtol/KK(1,1)); %5 times ransac threshold
%                    tmp = tmp & sqrt(sum((SIFT{i}.locs(index(1,:),[1,2])-SIFT{j}.locs(:,[1,2]))'.^2))'<settings.distance;
%                end
                
                matches(tmp) = index(1,tmp);
                fprintf('Found %d mathches\n',sum(matches ~= 0));
                
                ind2 = find(matches ~= 0);
                ind1 = matches(ind2);
                
                %Remove non uniqe matches (Caused by NN approximation?)             
                [~,I,~]=unique(ind1);
                ind2 = ind2(I);
                ind1 = ind1(I);
                
                pairwiseEst{i,j}.ind1 = ind1;
                pairwiseEst{i,j}.ind2 = ind2;
            else
                pairwiseEst{i,j}.ind1 = [];
                pairwiseEst{i,j}.ind2 = [];
            end
        else %No camera graph. Match all pairs.
            [i j length(seq)]
            
            
            tmp2 = normr(double(SIFT{j}.desc))';
            [index,dist] = vl_kdtreequery(kdtree,tmp1,tmp2,'MAXCOMPARISONS',50,'NUMNEIGHBORS',2);
            maxvals = sum(tmp1(:,index(1,:)).*tmp2,1)';
            secondmaxvals = sum(tmp1(:,index(2,:)).*tmp2,1)';
%                maxvals = sum(SIFT{i}.desc(index(1,:),:).*SIFT{j}.desc,2);
%                secondmaxvals = sum(SIFT{i}.desc(index(2,:),:).*SIFT{j}.desc,2);
                
            matches = zeros(1,size(SIFT{j}.desc,1));
                
            tmp = acos(maxvals) < distRatio*acos(secondmaxvals);
%            if ~isempty(settings.expectedF);
%                    F = settings.expectedF;
%                    p1 = pextend(normalize(SIFT{i}.locs(index(1,:),[2,1])',fc,cc,kc,alpha_c));
%                    p2 = pextend(normalize(SIFT{j}.locs(:,[2,1])',fc,cc,kc,alpha_c));
%                    l1 = p1'*F;l1norm=sqrt(l1(:,1).^2+l1(:,2).^2);d1 = abs(diag(l1 * p2)./l1norm);
%                    l2 = p2'*F';l2norm=sqrt(l2(:,1).^2+l2(:,2).^2);d2 = abs(diag(l2 * p1)./l2norm);
%                    dist=d1+d2;
%                    tmp = tmp & (dist < 5*settings.RANSAC_pixtol/KK(1,1)); %5 times ransac threshold
%                    tmp = tmp & sqrt(sum((SIFT{i}.locs(index(1,:),[1,2])-SIFT{j}.locs(:,[1,2]))'.^2))'<settings.distance;
%            end
               
            matches(tmp) = index(1,tmp);
            fprintf('Found %d mathches\n',sum(matches ~= 0));
            
            
            
            
%            [index,dist] = vl_kdtreequery(kdtree,SIFT{i}.desc',SIFT{j}.desc','MAXCOMPARISONS',40,'NUMNEIGHBORS',2);
%            maxvals = sum(SIFT{i}.desc(index(1,:),:).*SIFT{j}.desc,2);
%            secondmaxvals = sum(SIFT{i}.desc(index(2,:),:).*SIFT{j}.desc,2);
            
%            matches = zeros(1,size(SIFT{j}.desc,1));
%            matches(acos(maxvals) < distRatio*acos(secondmaxvals)) = index(1,acos(maxvals) < distRatio*acos(secondmaxvals));



            ind2 = find(matches ~= 0);
            ind1 = matches(ind2);
            
            %Remove non uniqe matches (Caused by NN approximation?)
            [~,I,~]=unique(ind1);
            ind2 = ind2(I);
            ind1 = ind1(I);

            pairwiseEst{i,j}.ind1 = ind1;
            pairwiseEst{i,j}.ind2 = ind2;
        end
     end
    end
end
%Matching done.
if settings.storesift == 0,
    %Only save the locations, not descriptors.
    for  i = 1:length(SIFT);
        SIFT{i}.desc=[];
    end
end

save(strcat(save_path,'pairwise_matchings.mat'), 'SIFT', 'pairwiseEst', 'imnames');

if 0 %Debug code
    i = 1;
    j = 2;
    
    if ~isempty(pairwiseEst{i,j})
        p1 = pextend(SIFT{i}.locs(pairwiseEst{i,j}.ind1,[2 1])');
        p2 = pextend(SIFT{j}.locs(pairwiseEst{i,j}.ind2,[2 1])');
        figure(1);
        plot(imagedata(strcat(img_path,imnames(i).name),p1),'numbered');
        figure(2);
        plot(imagedata(strcat(img_path,imnames(j).name),p2),'numbered');
    end    
end


if 0, %extract sift descriptors
    result = 'C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data/data1/result_data';
    load (result);
    pairwise_matching(settings); %put breakpoint after sift extraction and run above two lines again
    im_descriptors=[];
    for ii=1:length(u_uncalib.points),
        tmp = inv(settings.KK)*u_uncalib.points{ii};
        ss = uint8(zeros(size(tmp,2),128));
        for jj=1:size(tmp,2);
            fc = settings.KK([1 5]);
            cc = settings.KK(1:2,3);
            alpha_c = settings.KK(1,2)/fc(1);
            tmp2 = normalize(SIFT{ii}.locs(:,[2,1])',fc,cc,settings.kc,alpha_c);
            
            [res,resind]=min(abs(tmp2(2,:)-tmp(2,jj))+abs(tmp2(1,:)-tmp(1,jj)));
            if res>0.001,
                error('Not found');
            end
            ss(jj,:)=SIFT{ii}.desc(resind,:);
        end
        im_descriptors.sift{ii}=ss;
    end
    save(result,'P_uncalib','U','settings','u_uncalib','im_descriptors');
    
end