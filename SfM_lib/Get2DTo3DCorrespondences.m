% Takes as input a picture pic, the two matrices U and u saved in the
% red_str_mot.mat file by the sfm library, and the settings struct, and
% returns the 2D to 3D correspondences. The 2D correspondences are in
% normalized image coordinates. 

function [points2d, points3d] = Get2DTo3DCorrespondences(pic, U, u, settings)
    % Create matrices containing sift vectors and corresponding 3d points
    % for all 3d points. We do this because each 3d point will have several
    % descriptor vectors (one from each image it was seen in).  
    SIFT_vectors = zeros(128,size(U,2)); 
    count = zeros(1,size(U,2)); 
    for i = 1:length(u.sift)
        SIFT_vectors(:,u.index{i}) = SIFT_vectors(:,u.index{i}) + double(u.sift{i});
        count(:,u.index{i}) = count(:,u.index{i}) + 1; 
    end
    for i = 1:length(count)
        SIFT_vectors(:,i) = SIFT_vectors(:,i)/count(1,i); 
    end
    
    % Extract SIFT features from input picture
    pic = imresize(single(rgb2gray(pic)),settings.scale);
    [f,d] = vl_sift(pic, 'PeakThresh', settings.PeakThresh, 'EdgeThresh', settings.EdgeThresh);
    if (settings.rotationvariantSIFT==1)
        [f,d] = vl_sift(pic,'frames',[f(1:3,:);zeros(1,size(f,2))]);
    end
    locs = f([1,2],:)/settings.scale;
    
    % Match SIFT features in picture to those in 3D model
    [ind2, ind1] = matchsiftvectors(d, SIFT_vectors, settings.distRatio);
    
    % Return coordinates for matched image points in normalized coordinates
    fc = settings.KK([1 5]);
    cc = settings.KK(1:2,3);
    alpha_c = settings.KK(1,2)/fc(1);
    kc = settings.kc; 
    points2d = pextend(normalize(locs(:,ind2),fc,cc,kc,alpha_c));
    
    points3d = U(1:3,ind1);
end

