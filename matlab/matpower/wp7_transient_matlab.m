function wp7_transient_matlab()
%WP7_TRANSIENT_MATLAB  Transient stability and CCT, MATLAB side.
%
%   Run from the MATLAB app:
%       >> setup_paths
%       >> wp7_transient_matlab
%
%   Twin of python/studies/wp7_transient.py. Classical multi-machine swing
%   model, Kron-reduced to the internal machine nodes.
%
%   READ THIS BEFORE INTERPRETING ANY CCT FROM THIS FILE:
%
%   CCT is NOT comparable across scenarios with different numbers of dynamic
%   units. The stability criterion is rotor-angle separation, which only exists
%   for units that HAVE a rotor angle. Converting a machine to a grid-following
%   converter deletes it from the swing model, so the metric loses the failure
%   mode it exists to detect -- CCT then "improves" as grid-following
%   penetration rises, which is an artefact, not physics.
%
%   The controlled comparison is the virtual-inertia sweep: fixed topology,
%   fixed dispatch, fixed unit count, only the grid-forming inertia varies.
%   That is what part 3 does, and it is where the GFM claim is settled.
%
%   Writes results/tables/wp7_virtual_inertia_sweep_matlab.csv

define_constants;
cfg = gb_config();
outDir = fullfile(gb_root(), 'results', 'tables');
if ~exist(outDir, 'dir'), mkdir(outDir); end

scen = struct('case9', struct('fault', 7, 'trip', [6 7]), ...
              'case39', struct('fault', 16, 'trip', [16 17]));
names = fieldnames(scen);

fprintf('%s\n', repmat('=', 1, 78));
fprintf('WP7 -- transient stability and CCT (MATLAB)\n');
fprintf('%s\n', repmat('=', 1, 78));

%% ------------------------------------------------ 1. baseline validation
fprintf('\n1. Baseline CCT, all-synchronous system:\n');
fprintf('  %8s %7s %13s %7s %10s %10s\n', ...
        'case', 'fault', 'cleared by', 'units', 'H_total', 'CCT (s)');
for ci = 1:numel(names)
    cname = names{ci}; s = scen.(cname);
    mpc = ext2int(loadcase(cname));
    pf  = gb_newton(mpc);
    net = gb_dynamic_network(mpc, pf.V, [], []);
    kb  = findBranch(mpc, s.trip);
    r   = gb_cct(mpc, pf.V, net, s.fault, [], [], kb, 1.0);
    fprintf('  %8s %7d %8d-%-4d %7d %10.1f %10.3f\n', ...
        cname, s.fault, s.trip(1), s.trip(2), net.n, sum(net.H), r.cct);
end

%% ----------------------------------------------------- 2. the artefact
fprintf('\n2. Why a naive grid-following penetration sweep is INVALID:\n');
fprintf('  %8s %8s %9s %10s   verdict\n', 'pen%', 'rotors', 'H_total', 'CCT (s)');
cname = 'case39'; s = scen.(cname);
mpc = ext2int(loadcase(cname)); pf = gb_newton(mpc);
kb = findBranch(mpc, s.trip);
for pen = [0 40 60 80]
    [g, penActual] = gb_select_ibr(mpc, pen);
    net = gb_dynamic_network(mpc, pf.V, g, []);
    r = gb_cct(mpc, pf.V, net, s.fault, g, [], kb, 1.0);
    if pen == 0, verdict = 'valid baseline';
    else,        verdict = 'ARTEFACT -- do not report'; end
    fprintf('  %8.1f %8d %9.1f %10.3f   %s\n', ...
        penActual, net.n, sum(net.H), r.cct, verdict);
end
fprintf(['\n  CCT appears to improve. It has not: there are simply fewer rotors\n' ...
         '  left to separate. A grid-following converter destabilises through PLL\n' ...
         '  loss of synchronisation, which an electromechanical model cannot\n' ...
         '  represent -- an independent argument that EMT analysis is mandatory.\n']);

%% ------------------------------------------ 3. controlled inertia sweep
fprintf('\n3. CONTROLLED: fixed topology and unit count, only virtual inertia varies\n');
hValues = [0.5 1 2 4 8 16];
rows = struct('case',{},'pen_pct',{},'h_virtual_s',{},'h_total_s',{}, ...
              'n_dynamic_units',{},'cct_s',{});

for ci = 1:numel(names)
    cname = names{ci}; s = scen.(cname);
    mpc = ext2int(loadcase(cname)); pf = gb_newton(mpc);
    kb = findBranch(mpc, s.trip);
    seenPen = [];

    for pen = [60 80]
        [g, penActual] = gb_select_ibr(mpc, pen);
        if isempty(g), continue; end
        if any(abs(seenPen - penActual) < 0.05), continue; end
        seenPen(end+1) = penActual; %#ok<AGROW>
        gfmBuses = gb_split_gfm(mpc, g, 100);

        fprintf('\n  %s @ %.1f%% IBR, 100%% grid-forming\n', cname, penActual);
        fprintf('    %14s %10s %10s\n', 'H_virtual (s)', 'H_total', 'CCT (s)');
        first = NaN; last = NaN;
        for hv = hValues
            net = gb_dynamic_network(mpc, pf.V, g, gfmBuses, hv);
            r = gb_cct(mpc, pf.V, net, s.fault, g, gfmBuses, kb, 1.0);
            fprintf('    %14.1f %10.1f %10.3f\n', hv, sum(net.H), r.cct);
            if isnan(first), first = r.cct; end
            last = r.cct;
            rows(end+1) = struct('case', cname, 'pen_pct', penActual, ...
                'h_virtual_s', hv, 'h_total_s', sum(net.H), ...
                'n_dynamic_units', net.n, 'cct_s', r.cct); %#ok<AGROW>
        end
        if first > 0
            fprintf('    -> CCT %.3fs -> %.3fs (%.1fx)\n', first, last, last/first);
        end
    end
end

gb_writetable(rows, fullfile(outDir, 'wp7_virtual_inertia_sweep_matlab.csv'));
fprintf('\nWritten: %s\n', fullfile(outDir, 'wp7_virtual_inertia_sweep_matlab.csv'));
end


function k = findBranch(mpc, pair)
define_constants;
k = [];
for j = 1:size(mpc.branch, 1)
    if isequal(sort([mpc.branch(j, F_BUS) mpc.branch(j, T_BUS)]), sort(pair))
        k = j; return
    end
end
end
