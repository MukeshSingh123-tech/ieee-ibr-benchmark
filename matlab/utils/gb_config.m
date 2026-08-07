function cfg = gb_config(which)
%GB_CONFIG  Load the shared project configuration (cached).
%
%   cfg = GB_CONFIG()              % config/scenarios.yaml
%   cfg = GB_CONFIG('tolerances')  % config/tolerances.yaml
%
%   Reads the JSON exported by python/studies/export_config.py so that MATLAB
%   and Python are driven by the SAME scenario definitions. Never hard-code a
%   penetration level, fault location or tolerance in a .m file -- put it in
%   config/scenarios.yaml, re-export, and read it from here.

if nargin < 1 || isempty(which)
    which = 'scenarios';
end

persistent cache
if isempty(cache); cache = struct(); end
key = matlab.lang.makeValidName(which);
if isfield(cache, key)
    cfg = cache.(key);
    return
end

here = fileparts(mfilename('fullpath'));
f = fullfile(here, '..', '..', 'data', 'interchange', [which '.json']);
if exist(f, 'file') ~= 2
    error('gb_config:missing', ...
        ['Config bridge not found:\n  %s\n\n' ...
         'Generate it first:\n  cd <repo>/IEEE/python\n' ...
         '  python studies/export_config.py'], f);
end

cfg = jsondecode(fileread(f));
cache.(key) = cfg;
end
