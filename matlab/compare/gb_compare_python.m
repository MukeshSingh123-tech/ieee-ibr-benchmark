function gb_compare_python()
%GB_COMPARE_PYTHON  Diff the MATLAB results against the Python results.
%
%   Run from the MATLAB app AFTER both toolchains have produced tables:
%       >> setup_paths
%       >> wp1_baseline_matlab
%       >> wp3_faults_matlab
%       >> gb_compare_python
%
%   This is the step that makes the project a cross-tool benchmark rather than
%   two unrelated studies. Both sides read the SAME config/scenarios.yaml, solve
%   the SAME IEEE cases, and write the SAME column names -- so any disagreement
%   is a real difference in implementation, not in setup.
%
%   Pairs compared:
%       wp1_convergence.csv          <-> wp1_convergence_matlab.csv
%       wp1_crosstool.csv            <-> wp1_crosstool_matlab.csv
%       wp3_fault_error.csv          <-> wp3_fault_error_matlab.csv
%       wp5_gfm_wellposedness.csv    <-> wp5_gfm_wellposedness_matlab.csv
%
%   Writes results/tables/crosstool_comparison.csv

tableDir = fullfile(gb_root(), 'results', 'tables');

fprintf('%s\n', repmat('=', 1, 78));
fprintf('CROSS-TOOL COMPARISON -- Python/pandapower vs MATLAB/MATPOWER\n');
fprintf('%s\n', repmat('=', 1, 78));

rows = struct('study', {}, 'quantity', {}, 'n_compared', {}, ...
              'max_abs_diff', {}, 'mean_abs_diff', {}, 'tolerance', {}, 'verdict', {});

%% ---------------------------------------------------- WP1 iteration counts
rows = addComparison(rows, tableDir, ...
    'wp1_convergence.csv', 'wp1_convergence_matlab.csv', ...
    {'case', 'algorithm'}, 'iterations', 0, 'WP1 solver iterations', ...
    @(T) T(strcmp(string(T.source), 'gridbench'), :));

%% ------------------------------------------------ WP1 voltage agreement
rows = addComparison(rows, tableDir, ...
    'wp1_crosstool.csv', 'wp1_crosstool_matlab.csv', ...
    {'case', 'algorithm'}, 'max_dvm_pu', 1e-6, 'WP1 |dVm| vs reference solver', []);

%% ------------------------------------------------------- WP3 fault current
rows = addComparison(rows, tableDir, ...
    'wp3_fault_error.csv', 'wp3_fault_error_matlab.csv', ...
    {'case', 'bus', 'fault_type'}, 'if_classical_pu', 1e-6, ...
    'WP3 classical fault current (pu)', []);

rows = addComparison(rows, tableDir, ...
    'wp3_fault_error.csv', 'wp3_fault_error_matlab.csv', ...
    {'case', 'bus', 'fault_type'}, 'if_ibr_pu', 1e-3, ...
    'WP3 IBR-aware fault current (pu)', ...
    @(T) T(strcmp(string(T.converged), 'True'), :));

%% --------------------------------------------------- WP5 GFM mitigation
rows = addComparison(rows, tableDir, ...
    'wp5_gfm_wellposedness.csv', 'wp5_gfm_wellposedness_matlab.csv', ...
    {'case', 'gfm_share_pct'}, 'well_posed_pct', 1e-6, ...
    'WP5 well-posed fraction (%)', []);

%% ------------------------------------------------------------------ report
fprintf('\n%-38s %6s %13s %13s %11s %8s\n', ...
        'quantity', 'n', 'max |diff|', 'mean |diff|', 'tolerance', 'verdict');
for k = 1:numel(rows)
    r = rows(k);
    fprintf('%-38s %6d %13.3e %13.3e %11.1e %8s\n', ...
        r.quantity, r.n_compared, r.max_abs_diff, r.mean_abs_diff, ...
        r.tolerance, r.verdict);
end

nFail = sum(strcmp({rows.verdict}, 'FAIL'));
fprintf('\n%d of %d comparisons FAILED.\n', nFail, numel(rows));
if nFail > 0
    fprintf(2, ['A failure here means the two toolchains disagree on the same\n' ...
                'problem. Fix that before trusting any downstream result.\n']);
end

gb_writetable(rows, fullfile(tableDir, 'crosstool_comparison.csv'));
fprintf('\nWritten: %s\n', fullfile(tableDir, 'crosstool_comparison.csv'));
end


%% ======================================================================
function rows = addComparison(rows, dir_, pyFile, mlFile, keyCols, valueCol, ...
                              tol, label, filterFcn)
%ADDCOMPARISON  Join two result tables on keyCols and diff one value column.
pyPath = fullfile(dir_, pyFile);
mlPath = fullfile(dir_, mlFile);

if exist(pyPath, 'file') ~= 2 || exist(mlPath, 'file') ~= 2
    rows(end+1) = struct('study', pyFile, 'quantity', label, 'n_compared', 0, ...
        'max_abs_diff', NaN, 'mean_abs_diff', NaN, 'tolerance', tol, ...
        'verdict', 'SKIP');
    return
end

Tpy = readtable(pyPath, 'TextType', 'string');
Tml = readtable(mlPath, 'TextType', 'string');
if ~isempty(filterFcn)
    try, Tpy = filterFcn(Tpy); catch, end
    try, Tml = filterFcn(Tml); catch, end
end

if ~all(ismember([keyCols, {valueCol}], Tpy.Properties.VariableNames)) || ...
   ~all(ismember([keyCols, {valueCol}], Tml.Properties.VariableNames))
    rows(end+1) = struct('study', pyFile, 'quantity', label, 'n_compared', 0, ...
        'max_abs_diff', NaN, 'mean_abs_diff', NaN, 'tolerance', tol, ...
        'verdict', 'SKIP');
    return
end

kpy = joinKey(Tpy, keyCols);
kml = joinKey(Tml, keyCols);
[common, ipy, iml] = intersect(kpy, kml, 'stable');

if isempty(common)
    rows(end+1) = struct('study', pyFile, 'quantity', label, 'n_compared', 0, ...
        'max_abs_diff', NaN, 'mean_abs_diff', NaN, 'tolerance', tol, ...
        'verdict', 'SKIP');
    return
end

a = toNum(Tpy.(valueCol)(ipy));
b = toNum(Tml.(valueCol)(iml));
ok = isfinite(a) & isfinite(b);
d = abs(a(ok) - b(ok));

if isempty(d)
    verdict = 'SKIP'; mx = NaN; mn = NaN;
else
    mx = max(d); mn = mean(d);
    if mx <= tol, verdict = 'PASS'; else, verdict = 'FAIL'; end
end

rows(end+1) = struct('study', pyFile, 'quantity', label, ...
    'n_compared', sum(ok), 'max_abs_diff', mx, 'mean_abs_diff', mn, ...
    'tolerance', tol, 'verdict', verdict);
end


function k = joinKey(T, cols)
parts = strings(height(T), 1);
for i = 1:numel(cols)
    v = T.(cols{i});
    if isnumeric(v)
        parts = parts + "|" + string(round(v, 6));
    else
        parts = parts + "|" + string(v);
    end
end
k = parts;
end


function x = toNum(v)
if isnumeric(v), x = double(v(:)); else, x = str2double(string(v(:))); end
end
