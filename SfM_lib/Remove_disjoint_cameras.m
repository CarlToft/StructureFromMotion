function imnames = Remove_disjoint_cameras(settings)
imnames = settings.imnames;
save_path = settings.save_path;

load(strcat(save_path,'rotations_inliers.mat'));
load(strcat(save_path,'rotations.mat'));
load(strcat(save_path,'impoints.mat'));
load(strcat(save_path,'pairwise_geom.mat'));

estim_cameras = false(size(A));
for i = 1:length(A)
    estim_cameras(i) = ~isempty(A{i});
end

A = A(estim_cameras);
imnames = imnames(estim_cameras);
impoints.points = impoints.points(estim_cameras);
if settings.storesift==1,
    impoints.sift = impoints.sift(estim_cameras);
end
impoints.index = impoints.index(estim_cameras);
pairwise_geom = pairwise_geom(estim_cameras,estim_cameras);
maxinl = maxinl(estim_cameras,estim_cameras);

save(strcat(save_path,'rotations_inliers2.mat'),'maxinl');
save(strcat(save_path,'rotations2.mat'),'A');
save(strcat(save_path,'impoints2.mat'),'imnames','impoints');
save(strcat(save_path,'pairwise_geom2.mat'),'pairwise_geom');
