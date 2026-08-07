function modelName = build_ieee9_simscape(modelName, mode)
%BUILD_IEEE9_SIMSCAPE  Build the WSCC 9-bus network in Simscape Electrical.
%
%   build_ieee9_simscape()                 % discrete, 50 us  (EMT-style)
%   build_ieee9_simscape('gb_ieee9', 'phasor')
%
%   Run from the MATLAB app:
%       >> setup_paths
%       >> build_ieee9_simscape
%
%   Builds the ELECTRICAL network from the MATPOWER case9 data, so the Simulink
%   model and the phasor studies are driven by identical R, X, B and dispatch.
%   Requires Simscape Electrical (Specialized Power Systems).
%
%   WHY case9: it is the only test system that fits every licence cap in this
%   project (ETAP Student 15 buses, PSCAD Free 15 nodes), so the same network
%   can eventually be carried across all of them.
%
%   VALIDATION IS NOT OPTIONAL. After building, run to steady state and compare
%   the bus voltages against the MATPOWER load flow:
%
%       >> validate_ieee9_steady_state
%
%   If the pre-fault operating point does not match to ~1e-3 pu the network is
%   wired wrong, and every transient result taken from it is meaningless.

if nargin < 1 || isempty(modelName), modelName = 'gb_ieee9'; end
if nargin < 2 || isempty(mode),      mode = 'discrete';      end

if exist('powergui', 'file') ~= 4 && exist('powergui', 'file') ~= 2
    error('build_ieee9_simscape:noSPS', ...
        ['Simscape Electrical (Specialized Power Systems) not found.\n' ...
         'Check with:  license(''test'',''Simscape_Electrical'')\n' ...
         'Everything in WP1-WP8 runs without it; only the EMT work needs it.']);
end

cfg = gb_config();
Ts  = cfg.dynamics.emt_timestep_s;
f0  = cfg.meta.base_frequency_hz;

mpc = loadcase('case9');
mpc = ext2int(mpc);
define_constants;
baseMVA = mpc.baseMVA;

% reference load flow -- the model must reproduce this
pf = gb_newton(mpc);
if ~pf.converged
    error('build_ieee9_simscape:pf', 'reference load flow did not converge');
end

if bdIsLoaded(modelName), close_system(modelName, 0); end
new_system(modelName);
open_system(modelName);
set_param(modelName, 'StopTime', num2str(cfg.dynamics.sim_time_s), ...
                     'SolverType', 'Fixed-step', 'FixedStep', num2str(Ts));

%% --- powergui: REQUIRED, and must be added first ------------------------
add_block('powerlib/powergui', [modelName '/powergui'], 'Position', [30 30 110 70]);
if strcmpi(mode, 'phasor')
    set_param([modelName '/powergui'], 'SimulationMode', 'Phasor', ...
              'PhasorFrequency', num2str(f0));
else
    set_param([modelName '/powergui'], 'SimulationMode', 'Discrete', ...
              'SampleTime', num2str(Ts));
end

nb = size(mpc.bus, 1);
baseKV = mpc.bus(:, BASE_KV);
baseKV(baseKV <= 0) = 230;          % case9 transmission level

%% --- buses: a measurement block per bus ---------------------------------
x0 = 220; dy = 120;
for i = 1:nb
    name = sprintf('Bus%d', mpc.bus(i, BUS_I));
    blk = [modelName '/' name];
    add_block('powerlib/Measurements/Three-Phase V-I Measurement', blk, ...
        'Position', [x0, 40 + (i-1)*dy, x0 + 70, 100 + (i-1)*dy]);
    set_param(blk, 'VoltageMeasurement', 'phase-to-ground', ...
                   'UseLabels', 'on', ...
                   'VoltageLabel', sprintf('V%d', mpc.bus(i, BUS_I)), ...
                   'CurrentLabel', sprintf('I%d', mpc.bus(i, BUS_I)));
end

%% --- branches: PI sections from the MATPOWER R/X/B ----------------------
for k = 1:size(mpc.branch, 1)
    f = mpc.branch(k, F_BUS); t = mpc.branch(k, T_BUS);
    Zbase = (baseKV(f)^2) / baseMVA;              % ohms

    R = mpc.branch(k, BR_R) * Zbase;              % ohm
    L = mpc.branch(k, BR_X) * Zbase / (2*pi*f0);  % henry
    Bsh = mpc.branch(k, BR_B) / Zbase;            % siemens (total)
    C = Bsh / (2*pi*f0);                          % farad (total)

    name = sprintf('Line_%d_%d', f, t);
    blk = [modelName '/' name];
    add_block('powerlib/Elements/Three-Phase PI Section Line', blk, ...
        'Position', [x0 + 180, 40 + (k-1)*90, x0 + 280, 100 + (k-1)*90]);
    set_param(blk, 'Frequency', num2str(f0), ...
        'Resistances',  sprintf('[%g %g]', R, 3*R), ...
        'Inductances',  sprintf('[%g %g]', L, 3*L), ...
        'Capacitances', sprintf('[%g %g]', max(C,1e-12), max(C/2,1e-12)), ...
        'Length', '1');
end

%% --- loads: constant-impedance three-phase RLC --------------------------
for i = 1:nb
    Pd = mpc.bus(i, PD); Qd = mpc.bus(i, QD);
    if Pd == 0 && Qd == 0, continue; end
    name = sprintf('Load%d', mpc.bus(i, BUS_I));
    blk = [modelName '/' name];
    add_block('powerlib/Elements/Three-Phase Parallel RLC Load', blk, ...
        'Position', [x0 + 380, 40 + (i-1)*dy, x0 + 460, 100 + (i-1)*dy]);
    set_param(blk, 'Configuration', 'Y (grounded)', ...
        'NominalVoltage', num2str(baseKV(i)*1e3), ...
        'NominalFrequency', num2str(f0), ...
        'ActivePower', num2str(Pd*1e6), ...
        'InductivePower', num2str(max(Qd,0)*1e6), ...
        'CapacitivePower', num2str(max(-Qd,0)*1e6), ...
        'MeasurementLoad', 'off');
end

%% --- sources: synchronous machines at the generator buses ---------------
gbus = mpc.gen(:, GEN_BUS);
for g = 1:size(mpc.gen, 1)
    if mpc.gen(g, GEN_STATUS) <= 0, continue; end
    b = gbus(g);
    name = sprintf('Src%d', b);
    blk = [modelName '/' name];
    add_block('powerlib/Electrical Sources/Three-Phase Source', blk, ...
        'Position', [x0 - 160, 40 + (g-1)*dy, x0 - 80, 100 + (g-1)*dy]);

    i = find(mpc.bus(:, BUS_I) == b, 1);
    % Thevenin impedance from the SAME assumed short-circuit level the phasor
    % fault studies use, so the two toolchains agree on source strength
    scMVA = cfg.machine_data.slack.sc_multiple_of_load * sum(mpc.bus(:, PD));
    set_param(blk, 'Frequency', num2str(f0), ...
        'Vrms', num2str(pf.Vm(i) * baseKV(i) * 1e3), ...
        'Phase', num2str(pf.Va_deg(i)), ...
        'SpecifyImpedance', 'on', ...
        'BaseVoltage', num2str(baseKV(i)*1e3), ...
        'ShortCircuitLevel', num2str(scMVA*1e6), ...
        'XoverR', num2str(cfg.machine_data.slack.x_over_r));
end

%% --- fault block --------------------------------------------------------
fb = cfg.faults.locations.case9(1);
blk = [modelName '/Fault'];
add_block('powerlib/Elements/Three-Phase Fault', blk, ...
    'Position', [x0 + 560, 40, x0 + 640, 110]);
set_param(blk, 'FaultA', 'on', 'FaultB', 'on', 'FaultC', 'on', ...
    'GroundFault', 'on', ...
    'SwitchTimes', sprintf('[%g %g]', cfg.dynamics.fault_apply_s, ...
                            cfg.dynamics.fault_apply_s + 0.10), ...
    'FaultResistance', '1e-3', 'GroundResistance', '1e-3', ...
    'SnubberResistance', '1e6');

annotate(modelName, sprintf([ ...
    'WSCC 9-BUS -- built from MATPOWER case9 (%s mode, Ts = %g us)\n\n' ...
    'Wiring is left to you: connect Src -> Bus -> Line -> Bus -> Load, and\n' ...
    'the Fault block at bus %d. Block PARAMETERS are already set from the\n' ...
    'case data, which is the part that is tedious and easy to get wrong.\n\n' ...
    'THEN VALIDATE:  validate_ieee9_steady_state\n' ...
    'The pre-fault bus voltages must match the MATPOWER load flow to ~1e-3 pu.\n' ...
    'If they do not, the network is wired wrong and no transient result from\n' ...
    'it means anything.\n\n' ...
    'Reference load flow (MATPOWER):\n%s'], ...
    upper(mode), Ts*1e6, fb, referenceTable(mpc, pf)), [40 900]);

Simulink.BlockDiagram.arrangeSystem(modelName);
save_system(modelName, fullfile(gb_root(), 'matlab', 'simulink', [modelName '.slx']));
fprintf('Built %s (%s mode).\n', modelName, mode);
fprintf('Blocks are parameterised from case9; wire them, then run\n');
fprintf('validate_ieee9_steady_state before trusting any dynamic result.\n');
end


function s = referenceTable(mpc, pf)
define_constants;
lines = "";
for i = 1:size(mpc.bus, 1)
    lines = lines + sprintf('  Bus %d:  |V| = %.4f pu   angle = %+.2f deg\n', ...
        mpc.bus(i, BUS_I), pf.Vm(i), pf.Va_deg(i));
end
s = char(lines);
end


function annotate(m, text, pos)
a = Simulink.Annotation([m '/notes']);
a.Text = text;
a.Position = pos;
a.FontSize = 9;
end
