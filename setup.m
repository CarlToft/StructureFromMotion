% Reset path to avoid namespace clashes
restoredefaultpath;
lib_folder = fullfile('/','home','davidg','libs');

% Include local visionary folder
addpath(fullfile(pwd, 'visionary'));

% Include calibration toolbox folder
addpath(fullfile(lib_folder,'TOOLBOX_calib'));

% MOSEK
addpath(fullfile(lib_folder,'mosek','8','toolbox','r2014a'))

%Setup VLfeat.
fold = pwd;
cd (fullfile(lib_folder,'vlfeat', 'toolbox'))
vl_setup
cd(fold);

% CVX
cd (fullfile(lib_folder,'cvx'));
cvx_setup
cd (fold);