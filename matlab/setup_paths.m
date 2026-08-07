function setup_paths()
%SETUP_PATHS  Put the IEEE benchmark on the MATLAB path and check prerequisites.
%
%   Run this ONCE at the start of every MATLAB session:
%
%       >> cd('<...>/my-portfolio/IEEE/matlab')
%       >> setup_paths
%
%   It adds every subfolder of matlab/ to the path, then verifies that MATPOWER
%   is installed and that the JSON config exported from the Python side is
%   present. Both checks fail loudly with instructions rather than letting a
%   study die later with a confusing error.

here = fileparts(mfilename('fullpath'));
addpath(genpath(here));
fprintf('Added to path: %s (and subfolders)\n', here);

%% --- MATPOWER -----------------------------------------------------------
if exist('runpf', 'file') ~= 2
    fprintf(2, [ ...
        '\nMATPOWER NOT FOUND.\n' ...
        '  1. Download from https://matpower.org (MATPOWER 8.x)\n' ...
        '  2. Unzip somewhere permanent, e.g. C:\\matpower8\n' ...
        '  3. In MATLAB:  addpath(genpath(''C:\\matpower8''));  savepath\n' ...
        '  4. Verify:     test_matpower\n' ...
        '  5. Re-run setup_paths\n\n']);
else
    try
        v = mpver('all');
        fprintf('MATPOWER %s found.\n', v(1).Version);
    catch
        fprintf('MATPOWER found.\n');
    end
end

%% --- exported configuration ---------------------------------------------
cfgFile = fullfile(here, '..', 'data', 'interchange', 'scenarios.json');
if exist(cfgFile, 'file') ~= 2
    fprintf(2, [ ...
        '\nCONFIG NOT EXPORTED.\n' ...
        '  The MATLAB side reads the same scenarios.yaml the Python side does,\n' ...
        '  via a JSON bridge. Generate it once (and after any config change):\n\n' ...
        '      cd <repo>/IEEE/python\n' ...
        '      python studies/export_config.py\n\n']);
else
    fprintf('Config bridge present: %s\n', cfgFile);
end

%% --- output folders -----------------------------------------------------
for d = {'figures', 'tables', 'logs'}
    p = fullfile(here, '..', 'results', d{1});
    if ~exist(p, 'dir'); mkdir(p); end
end

fprintf('\nReady. Suggested first run:  wp1_baseline_matlab\n');
end
