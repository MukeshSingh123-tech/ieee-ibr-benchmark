function wp3_faults_matlab()
%WP3_FAULTS_MATLAB  Fault analysis and the grid-forming mitigation, MATLAB side.
%
%   Run from the MATLAB app:
%       >> setup_paths
%       >> wp3_faults_matlab
%
%   Mirrors python/studies/wp3_faults.py and wp5_mitigation.py so the two
%   toolchains can be diffed by gb_compare_python.
%
%   Three blocks:
%     1. Classical fault solve validated against the closed-form textbook
%        expressions -- if this does not match, nothing downstream is trustworthy.
%     2. Classical vs IBR-aware fault current across IBR penetration.
%     3. The grid-forming mitigation: does GFM restore a solvable problem?
%
%   Writes:
%       results/tables/wp3_fault_error_matlab.csv
%       results/tables/wp5_gfm_wellposedness_matlab.csv

define_constants;
cfg    = gb_config();
cases  = {'case14', 'case39'};
ftypes = {'3LG', 'SLG', 'LL', 'LLG'};
levels = cfg.penetration.levels_pct(:)';
outDir = fullfile(gb_root(), 'results', 'tables');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('%s\n', repmat('=', 1, 78));
fprintf('WP3 -- fault analysis with current-limited IBRs (MATLAB)\n');
fprintf('%s\n', repmat('=', 1, 78));

%% ---------------------------------------------- 1. closed-form validation
fprintf('\nClassical solver vs closed-form textbook expressions (IEEE 14-bus, bus 4):\n');
mpc = loadcase('case14');
pf  = gb_newton(mpc);
seq = gb_sequence_networks(mpc, pf.V);
k   = full(seq.e2i(4));
z1 = seq.Z1(k,k); z2 = seq.Z2(k,k); z0 = seq.Z0(k,k); vf = seq.Vprefault(k);

expected = struct('x3LG', abs(vf/z1), ...
                  'SLG',  abs(3*vf/(z1+z2+z0)), ...
                  'LL',   abs(sqrt(3)*vf/(z1+z2)));
names = {'3LG','SLG','LL'};
vals  = [expected.x3LG, expected.SLG, expected.LL];
for j = 1:3
    got = gb_fault(seq, 4, names{j}).iFaultMag;
    relErr = abs(got - vals(j)) / vals(j);
    fprintf('  %-4s closed form %10.6f   solver %10.6f   rel err %.2e  %s\n', ...
        names{j}, vals(j), got, relErr, ternary(relErr < 1e-9, 'PASS', 'FAIL'));
end

%% --------------------------------------- 2. classical vs IBR-aware sweep
rows = struct('case',{},'pen_pct',{},'bus',{},'fault_type',{}, ...
              'if_classical_pu',{},'if_ibr_pu',{},'error_pct',{}, ...
              'angle_shift_deg',{},'converged',{},'well_posed',{});

fprintf('\nFault current error, IBR-aware vs classical (converged solves only):\n');
for ci = 1:numel(cases)
    cname = cases{ci};
    mpc   = loadcase(cname);
    pf    = gb_newton(mpc);
    locs  = cfg.faults.locations.(cname)(:)';

    fprintf('\n  %s', cname);
    fprintf('%13s', ftypes{:}); fprintf('\n');

    for pen = levels
        [ibrGens, penActual] = gb_select_ibr(mpc, pen);
        if isempty(ibrGens), continue; end
        ibrBuses = unique(mpc.gen(ibrGens, GEN_BUS));
        seq = gb_sequence_networks(mpc, pf.V, ibrBuses);

        cells = zeros(1, numel(ftypes)); nOK = zeros(1, numel(ftypes));
        for fi = 1:numel(ftypes)
            errs = [];
            for b = locs
                a = gb_fault(seq, b, ftypes{fi});
                r = gb_fault(seq, b, ftypes{fi}, 'ibr_aware');
                if r.converged && a.iFaultMag > 1e-12
                    e = (r.iFaultMag - a.iFaultMag) / a.iFaultMag * 100;
                    errs(end+1) = e; %#ok<AGROW>
                end
                rows(end+1) = struct('case', cname, 'pen_pct', penActual, ...
                    'bus', b, 'fault_type', ftypes{fi}, ...
                    'if_classical_pu', a.iFaultMag, 'if_ibr_pu', r.iFaultMag, ...
                    'error_pct', (r.iFaultMag - a.iFaultMag)/a.iFaultMag*100, ...
                    'angle_shift_deg', gb_wrap_deg(r.iFaultAngleDeg - a.iFaultAngleDeg), ...
                    'converged', r.converged, 'well_posed', r.wellPosed); %#ok<AGROW>
            end
            if ~isempty(errs), cells(fi) = mean(errs); nOK(fi) = numel(errs); end
        end
        fprintf('  %5.1f%%', penActual);
        for fi = 1:numel(ftypes)
            if nOK(fi) > 0, fprintf('%12.2f%%', cells(fi));
            else,           fprintf('%13s', 'n/a'); end
        end
        fprintf('\n');
    end
end
gb_writetable(rows, fullfile(outDir, 'wp3_fault_error_matlab.csv'));

%% ------------------------------------- 3. grid-forming mitigation (WP5)
shares = [0 25 50 75 100];
rows2 = struct('case',{},'pen_pct',{},'gfm_share_pct',{}, ...
               'voltage_source_buses',{},'well_posed_pct',{},'converged_pct',{}, ...
               'mean_fault_current_pu',{});

fprintf('\n\nMITIGATION: does grid-forming control restore a solvable problem?\n');
for ci = 1:numel(cases)
    cname = cases{ci};
    mpc = loadcase(cname); pf = gb_newton(mpc);
    locs = cfg.faults.locations.(cname)(:)';

    for pen = [60 80 100]
        [ibrGens, penActual] = gb_select_ibr(mpc, pen);
        if isempty(ibrGens), continue; end
        ibrBuses = unique(mpc.gen(ibrGens, GEN_BUS));

        fprintf('\n  %s @ %.0f%% IBR penetration\n', cname, penActual);
        fprintf('    %10s %16s %12s %11s %14s\n', ...
                'GFM share', 'V-source buses', 'well-posed', 'converged', 'mean |If| pu');

        for s = shares
            gfmBuses = gb_split_gfm(mpc, ibrGens, s);
            seq = gb_sequence_networks(mpc, pf.V, ibrBuses, gfmBuses);

            nOK = 0; nWP = 0; nTot = 0; mags = [];
            for b = locs
                for fi = 1:numel(ftypes)
                    r = gb_fault(seq, b, ftypes{fi}, 'ibr_aware');
                    nTot = nTot + 1;
                    nOK = nOK + r.converged; nWP = nWP + r.wellPosed;
                    if r.converged, mags(end+1) = r.iFaultMag; end %#ok<AGROW>
                end
            end
            meanMag = NaN; if ~isempty(mags), meanMag = mean(mags); end
            fprintf('    %9d%% %16d %11.0f%% %10.0f%% %14.4f\n', ...
                s, numel(seq.voltageSourceBuses), nWP/nTot*100, nOK/nTot*100, meanMag);

            rows2(end+1) = struct('case', cname, 'pen_pct', penActual, ...
                'gfm_share_pct', s, ...
                'voltage_source_buses', numel(seq.voltageSourceBuses), ...
                'well_posed_pct', nWP/nTot*100, 'converged_pct', nOK/nTot*100, ...
                'mean_fault_current_pu', meanMag); %#ok<AGROW>
        end
    end
end
gb_writetable(rows2, fullfile(outDir, 'wp5_gfm_wellposedness_matlab.csv'));

fprintf('\nTables written to %s\n', outDir);
fprintf('Next: gb_compare_python\n');
end


function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
