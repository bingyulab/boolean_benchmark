%% Portable optPBN installation script
% This script installs optPBN and saves the path to user directory
% Works in any environment with proper $HOME paths

fprintf('Starting optPBN installation...\n');

% Use environment variable for portability
home_dir = getenv('HOME');
if isempty(home_dir)
    home_dir = '~';
end
optpbn_path = fullfile(home_dir, 'optPBN', 'optPBN_stand_alone_v2.2.3');
startup_dir = fullfile(home_dir, 'matlab_startup');

% Check if optPBN directory exists
if ~exist(optpbn_path, 'dir')
    error('optPBN directory not found: %s\nPlease install optPBN first!', optpbn_path);
end

% Create startup directory if it doesn't exist
if ~exist(startup_dir, 'dir')
    mkdir(startup_dir);
    fprintf('Created MATLAB startup directory: %s\n', startup_dir);
end

% Change to optPBN directory
cd(optpbn_path);

% Get the optPBN directory
fileposition = which('install');
optpbn_directory = regexprep(fileposition, 'install.m','');
fprintf('optPBN directory: %s\n', optpbn_directory);

%% 1) Add optimisation toolbox from SBToolbox2 to Matlab path
SBTB2_directory = [optpbn_directory filesep 'SBTOOLBOX2' filesep];
if exist(SBTB2_directory, 'dir')
    addpath(genpath(SBTB2_directory));
    fprintf('✓ Optimisation toolbox of SBToolbox2 is added to Matlab path\n');
else
    fprintf('⚠ Warning: SBTOOLBOX2 directory not found: %s\n', SBTB2_directory);
end

%% 2) Add BNPBN toolbox to Matlab path
BNPBN_directory = [optpbn_directory filesep 'pbn-matlab-toolbox' filesep];
if exist(BNPBN_directory, 'dir')
    addpath(BNPBN_directory);
    fprintf('✓ BNPBN toolbox is added to Matlab path\n');
else
    fprintf('⚠ Warning: BNPBN directory not found: %s\n', BNPBN_directory);
end

%% 3) Add optPBN toolbox to Matlab path
optPBN_directory = [optpbn_directory filesep 'optPBN_toolbox' filesep];
if exist(optPBN_directory, 'dir')
    addpath(optPBN_directory);
    fprintf('✓ optPBN scripts are added to Matlab path\n');
else
    fprintf('⚠ Warning: optPBN_toolbox directory not found: %s\n', optPBN_directory);
end

%% Save path to user directory
pathdef_file = fullfile(startup_dir, 'pathdef.m');
try
    savepath(pathdef_file);
    fprintf('✓ Path saved to user directory: %s\n', pathdef_file);
catch ME
    fprintf('⚠ Failed to save path: %s\n', ME.message);
end

%% Test installation
fprintf('\nTesting optPBN installation...\n');
optpbn_functions = {'add2estim', 'bestParams', 'evalCijAndSimulate'};
for i = 1:length(optpbn_functions)
    func_path = which(optpbn_functions{i});
    if ~isempty(func_path)
        fprintf('✓ %s: found\n', optpbn_functions{i});
    else
        fprintf('⚠ %s: not found\n', optpbn_functions{i});
    end
end
fprintf('\n✓ optPBN installation completed!\n');