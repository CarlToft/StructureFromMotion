function settings = reconstr_setup_corner

%Load calibration matrix (same calibration for all cameras at present)
load Calib_Results.mat KK
settings.KK = KK;
settings.kc = zeros(5,1);

%Select the images /home/toft/Documents/carltoft/tunnels/experiment3_data
settings.img_path = '/home/toft/Documents/carltoft/carl_olsson_system/corner/'; 
settings.imnames = dir(strcat(settings.img_path,'*.JPG'));

settings.expectedF = [];

settings.rotationvariantSIFT = 1;
settings.PeakThresh = 1;
settings.EdgeThresh = 10;
                       
                       
%path to where to save results
settings.save_path = './corner/';

%Path to Lowes SIFT-implementation
settings.SIFT_path = '/home/toft/Documents/carltoft/siftdemov4/';

%Setup VLfeat.
fold = pwd;
cd /home/toft/Documents/carltoft/vlfeat-0.9.13/toolbox/
vl_setup
cd(fold);

%Rescales the images to speed up SIFT.
settings.scale = 0.5;

% Same criterion as Lowe = 0.5 
settings.distRatio = 0.5; 

%RANSAC_threshold
settings.RANSAC_pixtol = 2.5; %Tolerans vid RANSAC-k�rning

%Minimum number matches to compute two-view geometris
settings.mininlnr = 20;

settings.storesift = 1; %save sift correspondences?

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

settings.forbidden={};
settings.epipoledistance = 150; %remove points around expected epipole