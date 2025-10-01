%% Portable MEIGO installation script
% This script installs MEIGO and saves the path to user directory
% Works in any environment with proper $HOME paths

fprintf('Starting MEIGO installation...\n');

% Use environment variable for portability
home_dir = getenv('HOME');
if isempty(home_dir)
    home_dir = '~';
end

meigo_path = fullfile(home_dir, 'MEIGO64', 'MEIGO');
startup_dir = fullfile(home_dir, 'matlab_startup');

% Check if MEIGO directory exists
if ~exist(meigo_path, 'dir')
    error('MEIGO directory not found: %s\nPlease install MEIGO first!', meigo_path);
end

% Create startup directory if it doesn't exist
if ~exist(startup_dir, 'dir')
    mkdir(startup_dir);
    fprintf('Created MATLAB startup directory: %s\n', startup_dir);
end

% Change to MEIGO directory
cd(meigo_path);

fprintf('MEIGO directory: %s\n', meigo_path);

%% Add MEIGO to Matlab path
addpath(genpath(meigo_path));
fprintf('✓ MEIGO is added to Matlab path\n');

%% Run MEIGO installer if it exists
try
    if exist('install_MEIGO.m', 'file')
        fprintf('Running MEIGO installer...\n');
        install_MEIGO;
        fprintf('✓ MEIGO installer completed\n');
    else
        fprintf('ℹ MEIGO installer not found, using manual setup\n');
    end
catch ME
    fprintf('⚠ MEIGO installer failed: %s\n', ME.message);
    fprintf('Continuing with manual path setup...\n');
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
fprintf('\nTesting MEIGO installation...\n');
meigo_functions = {'MEIGO', 'eSS'};
for i = 1:length(meigo_functions)
    func_path = which(meigo_functions{i});
    if ~isempty(func_path)
        fprintf('✓ %s: found at %s\n', meigo_functions{i}, func_path);
    else
        fprintf('⚠ %s: not found\n', meigo_functions{i});
    end
end

fprintf('\n✓ MEIGO installation completed!\n');
