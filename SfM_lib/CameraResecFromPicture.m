% Calculates camera pose given a picture (obtained using imread) and a 3d
% point cloud. The 3d cloud is given as the two inputs U and u that are
% saved in the red_str_mot.mat file created when running the sfm library.
% The settings struct used to create the 3d structure is needed as well.
% Resectioning is performed using 3D RANSAC. Should probably include bundle
% adjustment on camera pose as well. 
function P = CameraResecFromPicture(im, U, u, settings)
    [points2d, points3d] = Get2DTo3DCorrespondences(im, U, u, settings);
    points2d = normc(points2d); % convert image coords into unit vectors
    [P,inl] = CameraPoseRANSAC_Mex(points3d, points2d, 0.002)
end
