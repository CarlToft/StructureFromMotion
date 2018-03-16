function setup(lib_env);

% Reset path to avoid namespace clashes
restoredefaultpath;

if nargin<1,
    lib_env = getenv('SFM_LIB_PATH');
end

if isempty(lib_env)
    % No lib path was supplied.
    error('Please supply a path to the dependency libraries');
else
    lib_folder = fullfile(lib_env);
end

% Check that the path is absolute
is_relative = isempty(regexp(lib_folder(1),'[/\\]', 'once')) && ...
              isempty(regexp(lib_folder(2),':', 'once'));
if is_relative 
    lib_folder = fullfile(pwd, lib_folder);
    disp(['Correcting relative path to: ' lib_folder]);
end

% Include local folders
addpath(fullfile(pwd, 'visionary'));
addpath(fullfile(pwd, 'SfM_lib'));

mex('-outdir',...
    fullfile(pwd, 'visionary'),...
    fullfile(pwd, 'visionary','calibrated_fivepoint_helper.c'));

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