% numMatches = epipolarMatches(image1, image2, KK, kc)
% 
% Calculates the number of inliers for the epipolar geometry defined by the
% two input pictures. SIFT features are extracted from both images and
% matched, and then the number of inliers are found by performing 5pt
% RANSAC followed by bundle adjustment to calculate the essential matrix. 
% Both images are plotted together with lines connecting the found matches 
function numMatches = epipolarMatches(image1, image2, KK, kc)
    % Extract SIFT features 
    if (ndims(image1) == 3)
        image1 = single(rgb2gray(image1)); 
    else
        image1 = single(image1); 
    end
    %[f1,d1] = vl_sift(image1, 'PeakThresh', 1, 'EdgeThresh', 5);
    [f1,d1] = vl_dsift(image1, 'size', 8, 'step', 2);
    locs1 = f1([1,2],:); 
    
    if (ndims(image2) == 3)
        image2 = single(rgb2gray(image2)); 
    else
        image2 = single(image2);
    end
    %[f2,d2] = vl_sift(image2, 'PeakThresh', 1, 'EdgeThresh', 5);
    [f2,d2] = vl_dsift(image2, 'size', 8, 'step', 2);
    locs2 = f2([1,2],:); 
    
    % Match SIFT features between the two pictures 
    [ind1, ind2] = matchsiftvectors(d1,d2,0.8); 
    
    % Find matches in normalized image coordinates 
    fc = KK([1 5]);
    cc = KK(1:2,3);
    alpha_c = KK(1,2)/fc(1); 
    pp1 = pextend(normalize(locs1(:,ind1),fc,cc,kc,alpha_c));
    pp2 = pextend(normalize(locs2(:,ind2),fc,cc,kc,alpha_c)); 
    
    [inlierList,Pmax] = RANSAC_Essential(pp1, pp2, 0.004);
    numMatches = nnz(inlierList);
    
    %%%%%%%%%%% PLOT THE RESULTS AS WELL %%%%%%%%%%%% 
    imSize = size(image1);
    imbig = [image1, image2(1:imSize(1), 1:imSize(2))]; 
    figure; 
    imshow(uint8(imbig)); 
    im1Ind = ind1(find(inlierList));
    im2Ind = ind2(find(inlierList));
    for i = 1:length(im1Ind)
        line([locs1(1,im1Ind(i));locs2(1,im2Ind(i))+size(image1,2)], [locs1(2,im1Ind(i)); locs2(2,im2Ind(i))]);
        pause; 
    end
    
    figure; 
    plot(motion(Pmax)); 
end

