function lookups = create_lookups(U, u, cameras, nbr_pts_to_add)
    nbr_3d_pts = size(U, 2);
    nbr_cameras = size(u.points, 2);

    if nargin < 3
        % By default, create lookups for all cameras.
        cameras = 1:nbr_cameras;
    end
    if nargin < 4
        % By default, create lookups for all the points.
        % Indicated with -1.
        nbr_pts_to_add = -1;
    end

    lookups = struct();
    if nbr_pts_to_add == -1
        lookups.U_idx = 1 : nbr_3d_pts;
    else
        lookups.U_idx = nbr_3d_pts-nbr_pts_to_add+1 : nbr_3d_pts;
    end
    lookups.u_idx = cell(1, nbr_cameras);
    for j = cameras
        nbr_2d_pts = size(u.points{j}, 2);
        if nbr_pts_to_add == -1
            lookups.u_idx{j} = 1 : nbr_2d_pts;
        else
            lookups.u_idx{j} = nbr_2d_pts-nbr_pts_to_add+1 : nbr_2d_pts;
        end
    end
end
