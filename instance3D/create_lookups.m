function lookups = create_lookups(U, u, nbr_pts_to_add, start_idx)
    nbr_3d_pts = size(U, 2);
    nbr_views = size(u.points, 2);

    if nargin < 3
        % By default, create lookups for all the points.
        nbr_pts_to_add = nbr_3d_pts;
    end
    if nargin < 4
        % By default use the trailing points, i.e. start "nbr_pts_to_add" from the end.
        start_idx = nbr_3d_pts - nbr_pts_to_add + 1;
    end

    lookups = struct();
    lookups.U_idx = [];
    for j = start_idx : start_idx+nbr_pts_to_add-1
        lookups.U_idx = [lookups.U_idx, j];
    end
end
