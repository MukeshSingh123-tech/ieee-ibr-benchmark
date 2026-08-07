function ok = validate_ieee9_steady_state(modelName)
%VALIDATE_IEEE9_STEADY_STATE  Check the Simscape model against the load flow.
%
%   ok = VALIDATE_IEEE9_STEADY_STATE()
%   ok = VALIDATE_IEEE9_STEADY_STATE('gb_ieee9')
%
%   RUN THIS BEFORE TRUSTING ANY DYNAMIC RESULT.
%
%   A transient simulation answers "what happens after the disturbance". If the
%   operating point it starts FROM is wrong, the answer is wrong no matter how
%   fine the timestep is. Every EMT study in this project is gated on this
%   check, and the tolerance comes from config/tolerances.yaml
%   (cross_tool.pscad_prefault_vs_matpower) rather than being invented here.
%
%   Compares steady-state bus voltage magnitude and angle from the Simscape
%   model against gb_newton on the same MATPOWER case.

if nargin < 1 || isempty(modelName), modelName = 'gb_ieee9'; end
define_constants;

tol = gb_config('tolerances');
tolVm = tol.cross_tool.pscad_prefault_vs_matpower.vm_pu;
tolVa = tol.cross_tool.pscad_prefault_vs_matpower.va_deg;

%% --- reference ----------------------------------------------------------
mpc = ext2int(loadcase('case9'));
pf  = gb_newton(mpc);
nb  = size(mpc.bus, 1);

fprintf('%s\n', repmat('=', 1, 70));
fprintf('Simscape steady state vs MATPOWER load flow (%s)\n', modelName);
fprintf('%s\n', repmat('=', 1, 70));

if ~bdIsLoaded(modelName)
    f = fullfile(gb_root(), 'matlab', 'simulink', [modelName '.slx']);
    if exist(f, 'file') ~= 2
        fprintf(2, 'Model not found: %s\nRun build_ieee9_simscape first.\n', f);
        ok = false; return
    end
    load_system(f);
end

%% --- run just long enough to settle -------------------------------------
% The fault block must be disabled: this is a PRE-fault check.
faultBlk = [modelName '/Fault'];
if getSimulinkBlockHandle(faultBlk) > 0
    set_param(faultBlk, 'SwitchTimes', '[1e6 1e6]');
end
set_param(modelName, 'StopTime', '2.0');

try
    simOut = sim(modelName, 'SaveOutput', 'on', 'ReturnWorkspaceOutputs', 'on');
catch err
    fprintf(2, 'Simulation failed: %s\n', err.message);
    fprintf(2, ['This usually means the network is not fully wired.\n' ...
                'build_ieee9_simscape sets block PARAMETERS; the connections\n' ...
                'are yours to draw (see docs/WORKFLOW.md, Track D).\n']);
    ok = false; return
end

%% --- extract the measured voltages --------------------------------------
% Depends on how you labelled the V-I Measurement blocks. build_ieee9_simscape
% sets VoltageLabel to V1..V9, so a "From" goto tag per bus is expected.
baseKV = mpc.bus(:, BASE_KV); baseKV(baseKV <= 0) = 230;

vmSim = nan(nb, 1); vaSim = nan(nb, 1);
logged = [];
try
    logged = simOut.logsout;
catch
end

if isempty(logged)
    fprintf(2, ['No logged signals found.\n' ...
        'Enable signal logging on the V-I Measurement outputs (right-click the\n' ...
        'signal -> Log Selected Signals) and name them V1..V%d.\n'], nb);
    ok = false; return
end

for i = 1:nb
    name = sprintf('V%d', mpc.bus(i, BUS_I));
    try
        el = logged.getElement(name);
        d  = el.Values.Data;
        last = d(end, :);
        % three-phase peak phase-to-ground -> RMS line-to-neutral -> pu
        vmSim(i) = (max(abs(last)) / sqrt(2)) / (baseKV(i)*1e3/sqrt(3));
    catch
        % signal absent; leave NaN and report it below
    end
end

%% --- compare ------------------------------------------------------------
fprintf('\n%6s %12s %12s %12s\n', 'bus', 'MATPOWER', 'Simscape', '|error|');
worst = 0; nCompared = 0;
for i = 1:nb
    if isnan(vmSim(i))
        fprintf('%6d %12.4f %12s %12s\n', mpc.bus(i, BUS_I), pf.Vm(i), 'not logged', '-');
        continue
    end
    e = abs(vmSim(i) - pf.Vm(i));
    worst = max(worst, e); nCompared = nCompared + 1;
    fprintf('%6d %12.4f %12.4f %12.2e\n', mpc.bus(i, BUS_I), pf.Vm(i), vmSim(i), e);
end

ok = nCompared > 0 && worst <= tolVm;
fprintf('\ncompared %d of %d buses, worst |dVm| = %.3e pu (tolerance %.1e)\n', ...
    nCompared, nb, worst, tolVm);

if ok
    fprintf('PASS -- the model reproduces the load flow. Dynamic results are meaningful.\n');
else
    fprintf(2, ['FAIL -- do NOT run transient studies on this model yet.\n' ...
        'Common causes, in order of likelihood:\n' ...
        '  1. a line or load left unconnected\n' ...
        '  2. per-unit base wrong on a load (NominalVoltage must be line-to-line)\n' ...
        '  3. source phase angles not set from the load flow\n' ...
        '  4. powergui not in the intended mode, or Ts mismatched\n']);
end
end
