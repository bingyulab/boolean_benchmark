%% Complete Installation Script for optPBN and MEIGO
% This script installs both optPBN and MEIGO in any environment
% Prerequisites: optPBN and MEIGO64 directories must exist in $HOME

fprintf('=== Starting optPBN and MEIGO Installation ===\n\n');

% Use environment variable for portability
home_dir = getenv('HOME');
if isempty(home_dir)
    home_dir = '~';
end

startup_dir = fullfile(home_dir, 'matlab_startup');

% Create startup directory if it doesn't exist
if ~exist(startup_dir, 'dir')
    mkdir(startup_dir);
    fprintf('Created MATLAB startup directory: %s\n', startup_dir);
end

%% Install optPBN
fprintf('\n--- Installing optPBN ---\n');
optpbn_path = fullfile(home_dir, 'optPBN', 'optPBN_stand_alone_v2.2.3');

if exist(optpbn_path, 'dir')
    cd(optpbn_path);
    
    % Get the optPBN directory
    fileposition = which('install');
    optpbn_directory = regexprep(fileposition, 'install.m','');
    
    % Add optPBN components
    components = {
        {'SBTOOLBOX2', true};  % {folder, use genpath}
        {'pbn-matlab-toolbox', false};
        {'optPBN_toolbox', false}
    };
    
    for i = 1:length(components)
        comp_dir = [optpbn_directory filesep components{i}{1} filesep];
        if exist(comp_dir, 'dir')
            if components{i}{2}
                addpath(genpath(comp_dir));
            else
                addpath(comp_dir);
            end
            fprintf('✓ %s added to path\n', components{i}{1});
        else
            fprintf('⚠ %s directory not found\n', components{i}{1});
        end
    end
else
    fprintf('⚠ optPBN not found at: %s\n', optpbn_path);
end

%% Install MEIGO
fprintf('\n--- Installing MEIGO ---\n');
meigo_path = fullfile(home_dir, 'MEIGO64', 'MEIGO');

if exist(meigo_path, 'dir')
    cd(meigo_path);
    addpath(genpath(meigo_path));
    fprintf('✓ MEIGO added to path\n');
    
    % Run MEIGO installer if it exists
    try
        if exist('install_MEIGO.m', 'file')
            install_MEIGO;
            fprintf('✓ MEIGO installer completed\n');
        end
    catch ME
        fprintf('⚠ MEIGO installer warning: %s\n', ME.message);
    end
else
    fprintf('⚠ MEIGO not found at: %s\n', meigo_path);
end

%% Save paths
fprintf('\n--- Saving Paths ---\n');
pathdef_file = fullfile(startup_dir, 'pathdef.m');
try
    savepath(pathdef_file);
    fprintf('✓ Paths saved to: %s\n', pathdef_file);
catch ME
    fprintf('⚠ Failed to save paths: %s\n', ME.message);
end

%% Verification
fprintf('\n--- Verification ---\n');
test_functions = {
    {'optPBN', {'add2estim', 'bestParams', 'evalCijAndSimulate'}};
    {'MEIGO', {'MEIGO', 'eSS'}}
};

for pkg = 1:length(test_functions)
    pkg_name = test_functions{pkg}{1};
    functions = test_functions{pkg}{2};
    
    fprintf('%s functions:\n', pkg_name);
    found_count = 0;
    for i = 1:length(functions)
        func_path = which(functions{i});
        if ~isempty(func_path)
            fprintf('  ✓ %s\n', functions{i});
            found_count = found_count + 1;
        else
            fprintf('  ⚠ %s: not found\n', functions{i});
        end
    end
    
    if found_count == length(functions)
        fprintf('✅ %s: Fully installed (%d/%d functions)\n', pkg_name, found_count, length(functions));
    elseif found_count > 0
        fprintf('⚠ %s: Partially installed (%d/%d functions)\n', pkg_name, found_count, length(functions));
    else
        fprintf('❌ %s: Not installed\n', pkg_name);
    end
end

fprintf('\n🎉 Installation process completed!\n');
fprintf('💡 To use in future MATLAB sessions, run: startup\n');
fprintf('📁 Startup script location: %s/startup.m\n', startup_dir);
