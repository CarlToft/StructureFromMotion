% Reset path to avoid namespace clashes
restoredefaultpath;

lib_env = getenv('SFM_LIB_PATH');
if isempty(lib_env)
    % If you are NOT using the environment variable, 
    % change this match your folder structure.
    lib_folder = fullfile('/','home','davidg','libs');
else
    lib_folder = fullfile(lib_env);
end

% Include local visionary folder
addpath(fullfile(pwd, 'visionary'));

% Include calibration toolbox folder
addpath(fullfile(lib_folder,'TOOLBOX_calib'));

% MOSEK
addpath(fullfile(lib_folder,'mosek','8','toolbox','r2014a'))

% SIFT
addpath(fullfile(lib_folder,'siftDemoV4'));

%Setup VLfeat.
fold = pwd;
cd (fullfile(lib_folder,'vlfeat', 'toolbox'))
vl_setup
cd(fold);

% CVX
% cd (fullfile(lib_folder,'cvx'));
% cvx_setup
% cd (fold);

clear;