function modelName = build_frt_chart(modelName)
%BUILD_FRT_CHART  Create the IEEE 2800 fault-ride-through Stateflow chart.
%
%   modelName = BUILD_FRT_CHART()
%   modelName = BUILD_FRT_CHART('gb_frt')
%
%   Run from the MATLAB app:
%       >> setup_paths
%       >> build_frt_chart
%
%   Builds a Stateflow chart implementing the low-voltage ride-through logic
%   that IEEE Std 2800-2022 requires of a transmission-connected inverter. This
%   is genuinely discrete, event-driven behaviour -- states, entry actions,
%   guarded transitions and a timer -- which is exactly what Stateflow is for
%   and what a Simulink block diagram models badly.
%
%   States
%       NORMAL     V within band. Active power priority, normal current control.
%       LVRT       V below the ride-through threshold. The inverter MUST remain
%                  connected and switch to REACTIVE priority, injecting
%                  Iq = K1*(1 - V) to support the depressed voltage.
%       HVRT       V above the upper threshold. Absorb reactive power.
%       RECOVERY   Fault cleared. Active power ramps back at a limited rate --
%                  a step return would itself disturb the system.
%       TRIPPED    The ride-through envelope was exceeded. Disconnect and latch.
%
%   The envelope itself comes from config/scenarios.yaml (ibr.lvrt), so the
%   chart and the Python fault solver enforce the SAME curve.
%
%   Outputs of the chart drive the converter current references, so this chart
%   is what makes the Simulink converter model obey the standard rather than
%   just riding through everything unconditionally.

if nargin < 1 || isempty(modelName)
    modelName = 'gb_frt_chart';
end

cfg = gb_config();
lvrt = cfg.ibr.lvrt;                       % struct array: v_pu, t_s
iLimit = cfg.ibr.current_limit_pu;
k1 = 2.0;                                  % reactive support gain, pu/pu dip

%% ---------------------------------------------------------------- new model
if bdIsLoaded(modelName), close_system(modelName, 0); end
sfnew(modelName);
set_param(modelName, 'StopTime', num2str(cfg.dynamics.sim_time_s));

root  = sfroot;
chart = root.find('-isa', 'Stateflow.Chart', '-and', 'Path', modelName);
chart.Name = 'FRT';
chart.ActionLanguage = 'C';

%% ------------------------------------------------------------------- I/O
addData(chart, 'Vpu',       'Input',  'double');
addData(chart, 'Ipriority', 'Output', 'double');   % 0 = P priority, 1 = Q priority
addData(chart, 'Id_ref',    'Output', 'double');
addData(chart, 'Iq_ref',    'Output', 'double');
addData(chart, 'connected', 'Output', 'double');
addData(chart, 'state_id',  'Output', 'double');   % for logging/plots

addConst(chart, 'I_LIM',   iLimit);
addConst(chart, 'K1',      k1);
addConst(chart, 'V_LOW',   0.88);   % below this -> LVRT (IEEE 2800 continuous band)
addConst(chart, 'V_HIGH',  1.10);   % above this -> HVRT
addConst(chart, 'V_REC',   0.90);   % recovery threshold
addConst(chart, 'P_RAMP',  0.20);   % pu/s active power recovery ramp

% Ride-through envelope as a lookup: the maximum time permitted at each depth.
for k = 1:numel(lvrt)
    addConst(chart, sprintf('VRT_V%d', k), lvrt(k).v_pu);
    addConst(chart, sprintf('VRT_T%d', k), lvrt(k).t_s);
end
addConst(chart, 'N_VRT', numel(lvrt));

%% ----------------------------------------------------------------- states
sNormal = addState(chart, 'NORMAL', [30 30 200 90], sprintf([ ...
    'entry:\n' ...
    'Ipriority = 0;\n' ...            % active power priority
    'connected = 1;\n' ...
    'state_id = 0;\n' ...
    'during:\n' ...
    'Id_ref = I_LIM;\n' ...
    'Iq_ref = 0;']));

sLvrt = addState(chart, 'LVRT', [30 170 200 120], sprintf([ ...
    'entry:\n' ...
    'Ipriority = 1;\n' ...            % REACTIVE priority, per IEEE 2800
    'connected = 1;\n' ...
    'state_id = 1;\n' ...
    'during:\n' ...
    'Iq_ref = K1 * (1.0 - Vpu);\n' ...
    'if (Iq_ref > I_LIM) { Iq_ref = I_LIM; }\n' ...
    'Id_ref = sqrt(I_LIM*I_LIM - Iq_ref*Iq_ref);']));

sHvrt = addState(chart, 'HVRT', [300 170 200 120], sprintf([ ...
    'entry:\n' ...
    'Ipriority = 1;\n' ...
    'connected = 1;\n' ...
    'state_id = 2;\n' ...
    'during:\n' ...
    'Iq_ref = -K1 * (Vpu - 1.0);\n' ...
    'if (Iq_ref < -I_LIM) { Iq_ref = -I_LIM; }\n' ...
    'Id_ref = sqrt(I_LIM*I_LIM - Iq_ref*Iq_ref);']));

sRecov = addState(chart, 'RECOVERY', [300 30 200 90], sprintf([ ...
    'entry:\n' ...
    'Ipriority = 0;\n' ...
    'connected = 1;\n' ...
    'state_id = 3;\n' ...
    'during:\n' ...
    'Id_ref = Id_ref + P_RAMP * 0.001;\n' ...   % rate-limited restoration
    'if (Id_ref > I_LIM) { Id_ref = I_LIM; }\n' ...
    'Iq_ref = 0;']));

sTrip = addState(chart, 'TRIPPED', [165 340 200 80], sprintf([ ...
    'entry:\n' ...
    'connected = 0;\n' ...
    'Id_ref = 0;\n' ...
    'Iq_ref = 0;\n' ...
    'state_id = 4;']));

chart.defaultTransition = [];
addDefault(chart, sNormal);

%% ------------------------------------------------------------- transitions
addTrans(chart, sNormal, sLvrt,  '[Vpu < V_LOW]');
addTrans(chart, sNormal, sHvrt,  '[Vpu > V_HIGH]');
addTrans(chart, sLvrt,   sRecov, '[Vpu > V_REC]');
addTrans(chart, sHvrt,   sRecov, '[Vpu < V_HIGH]');
addTrans(chart, sRecov,  sNormal,'[Id_ref >= I_LIM - 0.001]');
addTrans(chart, sRecov,  sLvrt,  '[Vpu < V_LOW]');

% Ride-through envelope violation -> trip. `after(t, sec)` is the Stateflow
% temporal operator; the permitted duration depends on how deep the dip is, so
% one guarded transition per envelope point.
for k = 1:numel(lvrt)
    guard = sprintf('after(%g, sec)[Vpu < %g]', lvrt(k).t_s, lvrt(k).v_pu);
    addTrans(chart, sLvrt, sTrip, guard);
end

Simulink.BlockDiagram.arrangeSystem(modelName);
save_system(modelName, fullfile(gb_root(), 'matlab', 'stateflow', [modelName '.slx']));

fprintf('Built Stateflow FRT chart: %s\n', modelName);
fprintf('  states      : NORMAL, LVRT, HVRT, RECOVERY, TRIPPED\n');
fprintf('  envelope    : %d points from config/scenarios.yaml (ibr.lvrt)\n', numel(lvrt));
fprintf('  current lim : %.2f pu\n', iLimit);
fprintf('\nOpen it with:  open_system(''%s'')\n', modelName);
fprintf('Verify by feeding Vpu a step from 1.0 -> 0.2 -> 1.0 and checking\n');
fprintf('that state_id goes 0 -> 1 -> 3 -> 0 and never reaches 4 for a\n');
fprintf('120 ms dip (inside the 150 ms zero-voltage envelope).\n');
end


%% ======================================================================
%  Stateflow API helpers
%  ======================================================================
function d = addData(chart, name, scope, dtype)
d = Stateflow.Data(chart);
d.Name = name;
d.Scope = scope;
d.DataType = dtype;
end

function d = addConst(chart, name, value)
d = Stateflow.Data(chart);
d.Name = name;
d.Scope = 'Constant';
d.DataType = 'double';
d.Props.InitialValue = num2str(value, '%.10g');
end

function s = addState(chart, name, pos, label)
s = Stateflow.State(chart);
s.Name = name;
s.Position = pos;
s.LabelString = sprintf('%s\n%s', name, label);
end

function t = addTrans(chart, src, dst, label)
t = Stateflow.Transition(chart);
t.Source = src;
t.Destination = dst;
t.LabelString = label;
end

function t = addDefault(chart, dst)
t = Stateflow.Transition(chart);
t.Destination = dst;
t.SourceEndPoint = [dst.Position(1)+40, dst.Position(2)-30];
t.MidPoint       = [dst.Position(1)+40, dst.Position(2)-15];
end
