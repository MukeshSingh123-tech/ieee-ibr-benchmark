function modelName = build_current_limiter_chart(modelName)
%BUILD_CURRENT_LIMITER_CHART  Converter current-limit priority logic in Stateflow.
%
%   Run from the MATLAB app:
%       >> setup_paths
%       >> build_current_limiter_chart
%
%   A converter cannot exceed its current limit. When the commanded current
%   would, something has to give -- and WHICH component gives is a design
%   choice with direct consequences for grid stability:
%
%     ACTIVE PRIORITY    keep P, curtail Q. Maximises energy delivered. Correct
%                        in normal operation; wrong during a fault, because it
%                        withholds the reactive support that props the voltage up.
%
%     REACTIVE PRIORITY  keep Q, curtail P. Required during a fault by
%                        IEEE Std 2800-2022. This is what keeps a depressed
%                        voltage from collapsing further.
%
%     BALANCED           scale both by the same factor, preserving power factor.
%                        Simple, and wrong in both regimes.
%
%   The mode must switch at the moment the fault arrives and switch back
%   afterwards, which is discrete event-driven behaviour -- exactly what
%   Stateflow is for and what a block diagram of switches models badly.
%
%   Pair this with build_frt_chart: the FRT chart decides WHEN the converter is
%   in a fault, this chart decides WHAT IT DOES about the current limit.

if nargin < 1 || isempty(modelName), modelName = 'gb_current_limiter'; end

cfg = gb_config();
iLim = cfg.ibr.current_limit_pu;

if bdIsLoaded(modelName), close_system(modelName, 0); end
sfnew(modelName);
set_param(modelName, 'StopTime', num2str(cfg.dynamics.sim_time_s));

root  = sfroot;
chart = root.find('-isa', 'Stateflow.Chart', '-and', 'Path', modelName);
chart.Name = 'CurrentLimiter';
chart.ActionLanguage = 'C';

%% ------------------------------------------------------------------ I/O
addData(chart, 'Id_cmd',   'Input',  'double');   % commanded active current
addData(chart, 'Iq_cmd',   'Input',  'double');   % commanded reactive current
addData(chart, 'Vpu',      'Input',  'double');
addData(chart, 'frt_state','Input',  'double');   % from build_frt_chart
addData(chart, 'Id_out',   'Output', 'double');
addData(chart, 'Iq_out',   'Output', 'double');
addData(chart, 'limiting', 'Output', 'double');   % 1 when the limit binds
addData(chart, 'mode_id',  'Output', 'double');

addConst(chart, 'I_LIM', iLim);
addConst(chart, 'V_FAULT', 0.88);

%% --------------------------------------------------------------- states
%  Each state applies one priority rule. The magnitude check
%  sqrt(Id^2 + Iq^2) > I_LIM is what makes the converter NONLINEAR -- and it is
%  the same clamp the Python fault solver applies, so the two agree.
sActive = addState(chart, 'ACTIVE_PRIORITY', [40 40 260 150], sprintf([ ...
    'entry:\n' ...
    'mode_id = 0;\n' ...
    'during:\n' ...
    'mag = sqrt(Id_cmd*Id_cmd + Iq_cmd*Iq_cmd);\n' ...
    'if (mag > I_LIM) {\n' ...
    '  limiting = 1;\n' ...
    '  Id_out = Id_cmd;\n' ...
    '  if (Id_out > I_LIM) { Id_out = I_LIM; }\n' ...
    '  Iq_out = sqrt(I_LIM*I_LIM - Id_out*Id_out);\n' ...
    '  if (Iq_cmd < 0) { Iq_out = -Iq_out; }\n' ...
    '} else {\n' ...
    '  limiting = 0; Id_out = Id_cmd; Iq_out = Iq_cmd;\n' ...
    '}']));

sReactive = addState(chart, 'REACTIVE_PRIORITY', [340 40 260 150], sprintf([ ...
    'entry:\n' ...
    'mode_id = 1;\n' ...
    'during:\n' ...
    'mag = sqrt(Id_cmd*Id_cmd + Iq_cmd*Iq_cmd);\n' ...
    'if (mag > I_LIM) {\n' ...
    '  limiting = 1;\n' ...
    '  Iq_out = Iq_cmd;\n' ...
    '  if (Iq_out >  I_LIM) { Iq_out =  I_LIM; }\n' ...
    '  if (Iq_out < -I_LIM) { Iq_out = -I_LIM; }\n' ...
    '  Id_out = sqrt(I_LIM*I_LIM - Iq_out*Iq_out);\n' ...
    '} else {\n' ...
    '  limiting = 0; Id_out = Id_cmd; Iq_out = Iq_cmd;\n' ...
    '}']));

sBalanced = addState(chart, 'BALANCED', [190 250 260 130], sprintf([ ...
    'entry:\n' ...
    'mode_id = 2;\n' ...
    'during:\n' ...
    'mag = sqrt(Id_cmd*Id_cmd + Iq_cmd*Iq_cmd);\n' ...
    'if (mag > I_LIM) {\n' ...
    '  limiting = 1;\n' ...
    '  Id_out = Id_cmd * I_LIM / mag;\n' ...
    '  Iq_out = Iq_cmd * I_LIM / mag;\n' ...
    '} else {\n' ...
    '  limiting = 0; Id_out = Id_cmd; Iq_out = Iq_cmd;\n' ...
    '}']));

addLocal(chart, 'mag');
addDefault(chart, sActive);

%% ---------------------------------------------------------- transitions
% frt_state comes from build_frt_chart: 1 = LVRT, 2 = HVRT.
% The switch to reactive priority is the IEEE 2800 requirement.
addTrans(chart, sActive,   sReactive, '[frt_state == 1 || frt_state == 2]');
addTrans(chart, sReactive, sActive,   '[frt_state == 0]');

Simulink.BlockDiagram.arrangeSystem(modelName);
save_system(modelName, fullfile(gb_root(), 'matlab', 'stateflow', [modelName '.slx']));

fprintf('Built current-limiter chart: %s\n', modelName);
fprintf('  current limit : %.2f pu\n', iLim);
fprintf('  modes         : ACTIVE_PRIORITY, REACTIVE_PRIORITY, BALANCED\n');
fprintf('\nVerify: drive Id_cmd = 1.0, Iq_cmd = 0.8 (magnitude 1.28 > %.2f).\n', iLim);
fprintf('  In ACTIVE_PRIORITY   expect Id = 1.00, Iq = %.3f\n', sqrt(iLim^2 - 1.0));
fprintf('  In REACTIVE_PRIORITY expect Iq = 0.80, Id = %.3f\n', sqrt(iLim^2 - 0.64));
fprintf('  Set frt_state = 1 and confirm it switches.\n');
end


%% ======================================================================
function d = addData(chart, name, scope, dtype)
d = Stateflow.Data(chart);
d.Name = name; d.Scope = scope; d.DataType = dtype;
end

function d = addLocal(chart, name)
d = Stateflow.Data(chart);
d.Name = name; d.Scope = 'Local'; d.DataType = 'double';
end

function d = addConst(chart, name, value)
d = Stateflow.Data(chart);
d.Name = name; d.Scope = 'Constant'; d.DataType = 'double';
d.Props.InitialValue = num2str(value, '%.10g');
end

function s = addState(chart, name, pos, label)
s = Stateflow.State(chart);
s.Name = name; s.Position = pos;
s.LabelString = sprintf('%s\n%s', name, label);
end

function t = addTrans(chart, src, dst, label)
t = Stateflow.Transition(chart);
t.Source = src; t.Destination = dst; t.LabelString = label;
end

function t = addDefault(chart, dst)
t = Stateflow.Transition(chart);
t.Destination = dst;
t.SourceEndPoint = [dst.Position(1)+40, dst.Position(2)-30];
t.MidPoint       = [dst.Position(1)+40, dst.Position(2)-15];
end
