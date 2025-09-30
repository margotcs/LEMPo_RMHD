% setup.m
% 
% Add all necessary folders to the MATLAB path for this project.
% Run this script once at the start of your session, from any location.
%
% This setup will make sure that all source code and utilities are found,
% Usage:
%   In MATLAB, just type:
%       setup
%   or, if not in the repo root folder:
%       run('path/to/setup.m')
%

% Get the folder where setup.m is located
projectRoot = fileparts(mfilename('fullpath'));

% Add the main src folder (and all subfolders, e.g., utils)
srcPath = fullfile(projectRoot, 'src');
if exist(srcPath, 'dir')
    addpath(genpath(srcPath));
else
    warning('Could not find 00_src folder at expected location: %s', srcPath);
end


disp('MATLAB project setup complete. All relevant paths have been added.');