function settings = reconstr_setup

%Load calibration matrix (same calibration for all cameras at present)
load Calib_Results.mat KK;
settings.KK = KK;

settings.kc = zeros(5,1); %radiell distortion
settings.rotationvariantSIFT = 0;
settings.distRatio = 0.5; %Lowe's sift criterion
settings.expectedF = [];
settings.PeakThresh = 0;
settings.EdgeThresh = 10;
settings.forbidden = [];
settings.storesift = 0; %save sift correspondences?
settings.expectedEpipole = [];



%Select the images
settings.img_path = [pwd filesep];
settings.imnames = dir(fullfile(settings.img_path,'*.JPG'));

%path to where to save results
settings.save_path = '.\corner\';

%Path to Lowes SIFT-implementation
settings.SIFT_path = 'C:\Users\fredrik\Documents\MATLAB\siftdemov4\';

%Rescales the images to speed up SIFT.
settings.scale = 0.5;

%RANSAC_threshold
settings.RANSAC_pixtol = 2.5; %Tolerans vid RANSAC-körning

%Minimum number matches to compute two-view geometris
settings.mininlnr = 20;

%Minimum number of inliers to trust two-view results
settings.mincorrnr = 20;

%Threshold for inlier (in rotation averaging)
settings.roterrtol = 0.1; 

%Tolerans for first known-rotation run
settings.pixtol1 = 5; 

%Tolerans for second known-rotation run (don't know why I have two of
%these)
settings.pixtol2 = 5; 

%Points should be seen in at least this many cameras to be included in the
%optimization (to save memory).
settings.visviews = 3;

%Points in the final reconstruction that are uncertain in depth are removed
%after the reconstruction.
settings.uncertin_tol = settings.pixtol1; %Problems with the scale ambiguity?

%camera graph can be used to single out images to be matched if the
%sequence is not unordered.
settings.camera_graph = [];

%Has something to do with the point tracking. Sould probably alwasys be 1.
settings.merge_tracks = 1;
