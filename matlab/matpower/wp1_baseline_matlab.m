function wp1_baseline_matlab()
%WP1_BASELINE_MATLAB  WP1 on the MATLAB/MATPOWER side.
%
%   Run from the MATLAB app:
%       >> setup_paths
%       >> wp1_baseline_matlab
%
%   Does three things, mirroring python/studies/wp1_baseline.py exactly so the
%   two can be diffed by gb_compare_python:
%
%     1. Validates gb_ybus against MATPOWER's makeYbus (machine precision).
%     2. Runs MATPOWER's own runpf with NR / GS / FDXB / FDBX, and our
%        hand-written gb_newton / gb_gauss / gb_fdlf, on every test case.
%     3. Sweeps the Gauss-Seidel acceleration factor.
%
%   Writes:
%       results/tables/wp1_convergence_matlab.csv
%       results/tables/wp1_crosstool_matlab.csv
%       results/tables/wp1_gs_acceleration_matlab.csv

define_constants;
cases = {'case9', 'case14', 'case30', 'case39'};
outDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'tables');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('%s\n', repmat('=', 1, 78));
fprintf('WP1 -- classical load flow baseline (MATLAB / MATPOWER)\n');
fprintf('%s\n\n', repmat('=', 1, 78));

%% ---------------------------------------------------------------- Ybus gate
fprintf('Ybus validation (gb_ybus vs MATPOWER makeYbus):\n');
for c = cases
    mpc = loadcase(c{1});
    mpc = ext2int(mpc);
    Ymine = gb_ybus(mpc);
    Yref  = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
    err   = full(max(max(abs(Ymine - Yref))));
    fprintf('  %-8s max|dY| = %.3e', c{1}, err);
    if err < 1e-9, fprintf('   PASS\n'); else, fprintf('   FAIL\n'); end
end

%% -------------------------------------------------------------- convergence
algs = {'nr', 'gs', 'fdxb', 'fdbx'};
rows = struct('case', {}, 'n_bus', {}, 'algorithm', {}, 'source', {}, ...
              'converged', {}, 'iterations', {}, 'time_ms', {}, 'mismatch', {});

fprintf('\nSolver convergence -- hand-written (gb_*) vs MATPOWER runpf:\n');
fprintf('%8s %6s', 'case', 'n_bus');
for a = algs, fprintf('%12s', a{1}); end
fprintf('   (own / MATPOWER)\n');

for c = cases
    mpc = loadcase(c{1});
    nb  = size(mpc.bus, 1);
    line1 = sprintf('%8s %6d', c{1}, nb);

    for a = algs
        % --- our implementation ---
        switch a{1}
            case 'nr',   r = gb_newton(mpc);
            case 'gs',   r = gb_gauss(mpc);
            otherwise,   r = gb_fdlf(mpc, a{1});
        end
        rows(end+1) = struct('case', c{1}, 'n_bus', nb, 'algorithm', a{1}, ...
            'source', 'gridbench', 'converged', r.converged, ...
            'iterations', r.iterations, 'time_ms', r.elapsed_s*1e3, ...
            'mismatch', r.mismatch); %#ok<AGROW>

        % --- MATPOWER reference ---
        mpopt = mpoption('verbose', 0, 'out.all', 0, 'pf.alg', upper(a{1}));
        if strcmpi(a{1}, 'gs'), mpopt = mpoption(mpopt, 'pf.gs.max_it', 8000); end
        tic; rp = runpf(mpc, mpopt); tref = toc;
        rows(end+1) = struct('case', c{1}, 'n_bus', nb, 'algorithm', a{1}, ...
            'source', 'matpower', 'converged', rp.success, ...
            'iterations', NaN, 'time_ms', tref*1e3, 'mismatch', NaN); %#ok<AGROW>

        if r.converged
            line1 = [line1 sprintf('%12d', r.iterations)]; %#ok<AGROW>
        else
            line1 = [line1 sprintf('%12s', 'DIVERGED')];   %#ok<AGROW>
        end
    end
    fprintf('%s\n', line1);
end
gb_writetable(rows, fullfile(outDir, 'wp1_convergence_matlab.csv'));

%% ------------------------------------------------- hand-written vs MATPOWER
tol  = gb_config('tolerances');
rows2 = struct('case', {}, 'algorithm', {}, 'converged', {}, ...
               'max_dvm_pu', {}, 'max_dva_deg', {}, 'passes_gate', {});

fprintf('\nCross-check: gb_* solutions vs MATPOWER runpf (same case):\n');
for c = cases
    mpc   = loadcase(c{1});
    mpopt = mpoption('verbose', 0, 'out.all', 0);
    ref   = runpf(mpc, mpopt);
    refVm = ref.bus(:, VM);
    refVa = ref.bus(:, VA);

    for a = algs
        switch a{1}
            case 'nr',   r = gb_newton(mpc);
            case 'gs',   r = gb_gauss(mpc);
            otherwise,   r = gb_fdlf(mpc, a{1});
        end
        if strcmpi(a{1}, 'gs')
            gate = tol.cross_tool.handwritten_gs_vs_matpower;
        else
            gate = tol.cross_tool.handwritten_vs_matpower;
        end

        if r.converged
            dVm = max(abs(r.Vm - refVm));
            dVa = max(abs(gb_wrap_deg(r.Va_deg - refVa)));
            pass = dVm <= gate.vm_pu && dVa <= gate.va_deg;
        else
            dVm = NaN; dVa = NaN; pass = false;
        end
        rows2(end+1) = struct('case', c{1}, 'algorithm', a{1}, ...
            'converged', r.converged, 'max_dvm_pu', dVm, ...
            'max_dva_deg', dVa, 'passes_gate', pass); %#ok<AGROW>
        fprintf('  %-8s %-5s  |dVm|=%9.3e  |dVa|=%9.3e  %s\n', ...
            c{1}, a{1}, dVm, dVa, ternary(pass, 'PASS', 'FAIL'));
    end
end
gb_writetable(rows2, fullfile(outDir, 'wp1_crosstool_matlab.csv'));

%% ---------------------------------------------------- GS acceleration sweep
accels = [1.0 1.2 1.4 1.6 1.8];
rows3 = struct('case', {}, 'accel', {}, 'converged', {}, 'iterations', {}, ...
               'time_ms', {}, 'mismatch', {});

fprintf('\nGauss-Seidel iterations vs acceleration factor:\n');
fprintf('%8s', 'case'); fprintf('%10.1f', accels); fprintf('\n');
for c = cases
    mpc = loadcase(c{1});
    line = sprintf('%8s', c{1});
    for al = accels
        r = gb_gauss(mpc, struct('accel', al, 'maxIter', 8000));
        rows3(end+1) = struct('case', c{1}, 'accel', al, ...
            'converged', r.converged, 'iterations', r.iterations, ...
            'time_ms', r.elapsed_s*1e3, 'mismatch', r.mismatch); %#ok<AGROW>
        if r.converged
            line = [line sprintf('%10d', r.iterations)]; %#ok<AGROW>
        else
            line = [line sprintf('%10s', '---')];        %#ok<AGROW>
        end
    end
    fprintf('%s\n', line);
end
gb_writetable(rows3, fullfile(outDir, 'wp1_gs_acceleration_matlab.csv'));

fprintf('\nTables written to %s\n', outDir);
fprintf('Next: gb_compare_python  (diffs these against the Python results)\n');
end


function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
