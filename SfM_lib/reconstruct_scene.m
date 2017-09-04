function reconstruct_scene(datafolder)

thisfolder = pwd;

%Read settings
cd(datafolder)
settings = reconstr_setup;
cd(thisfolder)

%Perform pairwise matching for all pairs, (unless camera_graph is
%specified)
pairwise_matching(settings);

%Track points throughout the image collection. 
create_imdata(settings); 

%Compute all the relative orientations. Uses the 5 point algorithm for
%essential matrix.
rel_or_5points(settings);

%Performs rotation averaging using a RANSAC type procedure
%Does not perform very well for larger datasets due to the dimension of the
%search space. Will implement a new one.
settings.imnames = RANSAC_orientations(settings);

%Remove cameras that can't be estimated.
settings.imnames = Remove_disjoint_cameras(settings);

%Re-track points
create_imdata_2(settings); 

%Solve known rotation problem and run bundle.
settings.imnames = known_rotation_1(settings);

%Solve known rotation problem with updated rotations and run bundle.
settings.imnames = known_rotation_2(settings);

%Re-triangulates points with current estmated cameras
%Might increase number of inliers slightly.
triang_outl_reconst(settings);

%Remove points that are uncertain in the depth direction.
[U,P_uncalib,u_uncalib] = remove_uncertin_points(settings); 
save(fullfile(settings.save_path,'str_mot3.mat'), 'U', 'P_uncalib', 'u_uncalib');
save(fullfile(settings.save_path,'settings3.mat'), 'settings');

%Plot result
plot_result(settings,P_uncalib,U,u_uncalib);
