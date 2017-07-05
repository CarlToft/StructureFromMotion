% [U, inliers] = triangulate_pt_ransac(xlist, Plist, iters, threshold)
%
% Robust triangulation of single point seen in multiple cameras using
% RANSAC. The points are triangulated using midpoint triangulation from two
% views inside the RANSAC loop. 
% 
% INPUTS:           xlist:       2xN matrix of image coordinates in
%                                normalized coordinates
%                   Plist:       Nx1 cell array of normalized camera
%                                matrices
%                   iters:       number of RANSAC iterations 
%                   threshold:   maximum allowed reprojection error for
%                                inliers
% 
% OUTPUTS:          U:           [4x1] vector containing best triangulated
%                                3d point
%                   inliers:     [Nx1] vector of booleans denoting whether
%                                the returned U is an inlier in each of the
%                                cameras in Plist
% 
function [U, inliers] = triangulate_pt_ransac(xlist, Plist, iters, threshold)
    N = length(Plist); 
    
    inliers = zeros(N,1); 
    bestInliers = zeros(N,1);
    bestNumInliers = 0; 
    bestU = zeros(4,1); 
    
    xlist = pextend(xlist); 
    for ii = 1:iters
        cameras = randperm(N,2); 
        U = intsec2views_midpoint(Plist{cameras(1)}, Plist{cameras(2)}, xlist(:,cameras(1)), xlist(:,cameras(2))); 
        
        % Reproject the triangulated points in all other views 
        proj = zeros(3,N); 
        for k = 1:N
            proj(:,k) = pflat(Plist{k}*U);
            inliers(k) = (norm(xlist(:,k)-proj(:,k)) < threshold); 
        end
        
        if (nnz(inliers) > bestNumInliers)
            bestInliers = inliers; 
            bestNumInliers = nnz(bestInliers); 
            bestU = U; 
        end
    end
    
    U = bestU; 
    inliers = bestInliers; 
end

