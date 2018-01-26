function [P, inliers] = ransac3pts_weighted(u, U, weights, tol, iters)
    % [P, inliers] = ransac3pts_weighted(u, U, weights, tol, iters)
    %
    % Performs weighted 3-point-RANSAC for pose of calibrated camera. u
    % should be a 2xN matrix with normalized image coordinates, U should be
    % a 3xN matrix with corresponding 3D points. weights is a vector of
    % length N where element i denotes the probability of picking
    % correspondence i in the RANSAC loop. tol is the reprojection
    % threshold for correspondences to be considered inliers, measured in
    % the z=1 image plane. iters denotes the number of ransac iterations. 
    
    N = size(u,2); 
    U_h = pextend(U); 
    P = zeros(3,4); 
    inliers = zeros(N,1); 
    for iter = 1:iters
        % Choose three points, make sure the correspondences are are unique
        % (sometimes sift triggers twice at the same point), and fit camera
        % pose to correspondences
        perm = randsample(N,3,true,weights); 
        if (all(U(:,perm(1))-U(:,perm(2)) == 0) || all(U(:,perm(2))-U(:,perm(3)) == 0) || all(U(:,perm(1))-U(:,perm(3)) == 0))
            continue;
        end
        if (all(u(:,perm(1))-u(:,perm(2)) == 0) || all(u(:,perm(2))-u(:,perm(3)) == 0) || all(u(:,perm(1))-u(:,perm(3)) == 0))
            continue;
        end
        
        motlist = resec3pts(imagedata([],u(:,perm)),structure(U(:,perm)));
        P_ransac = getcameras(motlist);
        
        % Go through all found cameras and count inliers
        for k = 1:length(P_ransac)
            U_proj = P_ransac{k}*U_h; 
            inFront = (U_proj(3,:) > 0);
            U_proj = pflat(U_proj); 
            currInliers = inFront & (sum((U_proj(1:2,:)-u).^2,1) < tol*tol);
             
            if (nnz(currInliers) > nnz(inliers))
                inliers = currInliers;
                P = P_ransac{k};
            end
        end
        
        if (mod(iter, 250) == 0)
            % display(['Finished iteration ', num2str(iter), '!']);
        end
    end
end

