function modelName = build_converter_control(mode, modelName)
%BUILD_CONVERTER_CONTROL  Build a GFL or GFM converter control model in Simulink.
%
%   build_converter_control('gfl')        % grid-following
%   build_converter_control('gfm')        % grid-forming
%
%   Run from the MATLAB app:
%       >> setup_paths
%       >> build_converter_control('gfl')
%       >> build_converter_control('gfm')
%
%   Builds the CONTROL side of the converter from plain Simulink blocks (Sources,
%   Math, Continuous, Discontinuities). It is kept separate from the Simscape
%   power network on purpose: control logic built this way is portable and can be
%   unit-tested on its own, whereas the electrical network is far easier to wire
%   in the GUI (see docs/WORKFLOW_MATLAB.md, step 5).
%
%   GRID-FOLLOWING (GFL)
%       Synchronises to the grid with a PLL and injects a CURRENT.
%       It needs a stiff voltage to lock onto, so it degrades as SCR falls --
%       this is the model whose PLL goes unstable in weak grids, which is what
%       makes SCR >= 3 an interconnection screen.
%
%       theta   <- PLL(Vabc)
%       Id_ref  <- P droop / MPPT
%       Iq_ref  <- voltage support (from the FRT Stateflow chart)
%       Vdq_ref <- PI(Idq_ref - Idq) + feedforward
%
%   GRID-FORMING (GFM)
%       Behaves as a VOLTAGE source behind an impedance and sets its own angle
%       from a power droop -- no PLL, so no PLL instability. It can operate at
%       SCR < 1 and can black-start, which a GFL converter cannot.
%
%       omega <- omega0 + mp*(Pset - P)        active power / frequency droop
%       theta <- integral(omega)
%       E     <- E0    + mq*(Qset - Q)         reactive power / voltage droop
%       Vdq_ref <- E - Zvirtual * Idq          virtual impedance for current limiting

if nargin < 1 || isempty(mode),  mode = 'gfl'; end
if nargin < 2 || isempty(modelName)
    modelName = ['gb_' lower(mode) '_control'];
end
mode = lower(mode);

cfg = gb_config();
Ts  = cfg.dynamics.emt_timestep_s;

if bdIsLoaded(modelName), close_system(modelName, 0); end
new_system(modelName);
open_system(modelName);
set_param(modelName, 'StopTime', num2str(cfg.dynamics.sim_time_s), ...
                     'SolverType', 'Fixed-step', ...
                     'FixedStep',  num2str(Ts));

switch mode
    case 'gfl'
        buildGFL(modelName, cfg);
    case 'gfm'
        buildGFM(modelName, cfg);
    otherwise
        error('build_converter_control:mode', 'mode must be ''gfl'' or ''gfm''');
end

Simulink.BlockDiagram.arrangeSystem(modelName);
outFile = fullfile(gb_root(), 'matlab', 'simulink', [modelName '.slx']);
save_system(modelName, outFile);
fprintf('Built %s control model: %s\n', upper(mode), outFile);
end


%% ======================================================================
function buildGFL(m, cfg)
%BUILDGFL  Grid-following: PLL + dq current control.

pll = cfg.dynamics.gfl.pll_bandwidth_hz;
ci  = cfg.dynamics.gfl.current_loop_bandwidth_hz;
f0  = cfg.meta.base_frequency_hz;

% --- inputs ---
addBlk(m, 'simulink/Sources/In1',  'Vabc',    [40  40  70  60]);
addBlk(m, 'simulink/Sources/In1',  'Iabc',    [40 120  70 140]);
addBlk(m, 'simulink/Sources/In1',  'Id_ref',  [40 200  70 220]);
addBlk(m, 'simulink/Sources/In1',  'Iq_ref',  [40 260  70 280]);
set_param([m '/Iabc'],   'Port', '2');
set_param([m '/Id_ref'], 'Port', '3');
set_param([m '/Iq_ref'], 'Port', '4');

% --- PLL: a second-order loop whose bandwidth is the key weak-grid parameter.
% Higher bandwidth tracks faster but is exactly what destabilises at low SCR.
kp_pll = 2 * 2*pi*pll;
ki_pll = (2*pi*pll)^2;
addSub(m, 'PLL', [140 30 260 90], sprintf([ ...
    'PLL, %.0f Hz bandwidth.\\n' ...
    'Kp = %.4g, Ki = %.4g\\n' ...
    'Locks theta to Vabc. In a weak grid the PLL sees its OWN\\n' ...
    'current affect the terminal voltage, closing a positive\\n' ...
    'feedback loop -- the standard GFL weak-grid instability.'], pll, kp_pll, ki_pll));

addSub(m, 'abc_to_dq_V', [140 110 260 160], 'Park transform on Vabc using theta');
addSub(m, 'abc_to_dq_I', [140 180 260 230], 'Park transform on Iabc using theta');

% --- current controller: PI on Id, Iq with cross-coupling decoupling ---
kp_i = 2*pi*ci * 0.001;      % L = 1 mH assumed; retune with the real filter
ki_i = kp_i * 2*pi*ci / 10;
addSub(m, 'CurrentControl_dq', [320 120 470 230], sprintf([ ...
    'PI current control, %.0f Hz bandwidth\\n' ...
    'Kp = %.4g, Ki = %.4g\\n' ...
    'Vd_ref = PI(Id_ref-Id) - wL*Iq + Vd\\n' ...
    'Vq_ref = PI(Iq_ref-Iq) + wL*Id + Vq'], ci, kp_i, ki_i));

addSub(m, 'CurrentLimiter', [520 150 650 210], sprintf([ ...
    'Magnitude limit |I| <= %.2f pu with priority selection.\\n' ...
    'Priority comes from the FRT Stateflow chart:\\n' ...
    '  0 = active priority  (normal operation)\\n' ...
    '  1 = reactive priority (during a fault, per IEEE 2800)\\n' ...
    'This block is where the converter stops being linear.'], cfg.ibr.current_limit_pu));

addSub(m, 'dq_to_abc', [700 150 820 210], 'Inverse Park -> Vabc_ref for the bridge');

addBlk(m, 'simulink/Sinks/Out1', 'Vabc_ref', [870 160 900 180]);
addBlk(m, 'simulink/Sinks/Out1', 'theta',    [870 40  900 60]);
set_param([m '/theta'], 'Port', '2');

annotate(m, sprintf([ ...
    'GRID-FOLLOWING CONVERTER (GFL)\n' ...
    'Base frequency %g Hz, timestep %g us\n\n' ...
    'Injects current, synchronised by a PLL. Requires a stiff grid:\n' ...
    'stability degrades as SCR falls and is the reason SCR >= 3 is\n' ...
    'used as an interconnection screen. Contributes NO inertia.\n\n' ...
    'Connect CurrentLimiter/priority to the FRT chart output\n' ...
    '(build_frt_chart) to make this model IEEE 2800 compliant.'], ...
    f0, cfg.dynamics.emt_timestep_s*1e6), [40 340]);
end


%% ======================================================================
function buildGFM(m, cfg)
%BUILDGFM  Grid-forming: power droop + virtual impedance, no PLL.

mp = cfg.dynamics.gfm.p_droop_pct / 100;
mq = cfg.dynamics.gfm.q_droop_pct / 100;
zv = cfg.dynamics.gfm.virtual_impedance_pu;
H  = cfg.dynamics.gfm.inertia_const_s;
f0 = cfg.meta.base_frequency_hz;

addBlk(m, 'simulink/Sources/In1', 'Vabc', [40  40  70  60]);
addBlk(m, 'simulink/Sources/In1', 'Iabc', [40 110  70 130]);
addBlk(m, 'simulink/Sources/In1', 'Pset', [40 180  70 200]);
addBlk(m, 'simulink/Sources/In1', 'Qset', [40 240  70 260]);
set_param([m '/Iabc'], 'Port', '2');
set_param([m '/Pset'], 'Port', '3');
set_param([m '/Qset'], 'Port', '4');

addSub(m, 'PowerCalc', [140 60 260 130], 'Instantaneous P, Q from Vabc, Iabc (filtered)');

addSub(m, 'ActiveDroop', [320 40 470 110], sprintf([ ...
    'ACTIVE POWER / FREQUENCY DROOP\\n' ...
    'omega = omega0 + mp*(Pset - P),  mp = %.4g\\n' ...
    'theta = integral(omega)\\n\\n' ...
    'NO PLL. The converter SETS the angle instead of tracking it,\\n' ...
    'which is why a GFM unit is stable at SCR < 1 and can\\n' ...
    'black-start, and why GFL cannot.'], mp));

addSub(m, 'ReactiveDroop', [320 150 470 220], sprintf([ ...
    'REACTIVE POWER / VOLTAGE DROOP\\n' ...
    'E = E0 + mq*(Qset - Q),  mq = %.4g'], mq));

if strcmpi(cfg.dynamics.gfm.control, 'vsm')
    addSub(m, 'VirtualInertia', [320 250 470 320], sprintf([ ...
        'VIRTUAL SYNCHRONOUS MACHINE\\n' ...
        '2H d(omega)/dt = Pset - P - D*(omega-omega0),  H = %.2f s\\n\\n' ...
        'Emulates the swing equation, so the converter contributes\\n' ...
        'SYNTHETIC INERTIA and limits RoCoF -- addressing exactly the\\n' ...
        'problem quantified in the WP5 inertia study.'], H));
end

addSub(m, 'VirtualImpedance', [520 120 670 200], sprintf([ ...
    'VIRTUAL IMPEDANCE  Zv = %.3g + j%.3g pu\\n' ...
    'Vref = E - Zv*I\\n\\n' ...
    'Provides current limiting WITHOUT abandoning voltage-source\\n' ...
    'behaviour. A GFM unit that simply clamps its current reverts to\\n' ...
    'grid-following during the fault and loses its advantage.'], zv.r, zv.x));

addSub(m, 'dq_to_abc', [720 130 830 190], 'Inverse Park -> Vabc_ref');

addBlk(m, 'simulink/Sinks/Out1', 'Vabc_ref', [880 150 910 170]);
addBlk(m, 'simulink/Sinks/Out1', 'theta',    [880 60  910 80]);
set_param([m '/theta'], 'Port', '2');

annotate(m, sprintf([ ...
    'GRID-FORMING CONVERTER (GFM) -- %s control\n' ...
    'Base frequency %g Hz, timestep %g us\n\n' ...
    'Acts as a voltage source behind an impedance and sets its own\n' ...
    'angle from power droop. Stable in weak grids, can black-start,\n' ...
    'and (with VSM) supplies synthetic inertia.\n\n' ...
    'The WP12 study compares CCT for GFM vs GFL at matched penetration.'], ...
    upper(cfg.dynamics.gfm.control), f0, cfg.dynamics.emt_timestep_s*1e6), [40 360]);
end


%% ======================================================================
function h = addBlk(m, src, name, pos)
h = add_block(src, [m '/' name], 'Position', pos);
end

function h = addSub(m, name, pos, note)
%ADDSUB  Placeholder subsystem carrying its specification as an annotation.
%   The maths is documented in the block so the model is self-describing;
%   fill in the internals following docs/WORKFLOW_MATLAB.md.
h = add_block('simulink/Ports & Subsystems/Subsystem', [m '/' name], 'Position', pos);
delete_line([m '/' name], 'In1/1', 'Out1/1');
delete_block([m '/' name '/In1']);
delete_block([m '/' name '/Out1']);
add_block('built-in/Note', [m '/' name '/spec'], 'Position', [40 40], ...
          'Text', note);
end

function annotate(m, text, pos)
a = Simulink.Annotation([m '/note']);
a.Text = text;
a.Position = pos;
a.FontSize = 10;
end
