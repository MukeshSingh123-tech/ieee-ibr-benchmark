function modelName = build_relay_chart(modelName)
%BUILD_RELAY_CHART  Protection relay logic (pickup / time / trip / reclose).
%
%   Run from the MATLAB app:
%       >> setup_paths
%       >> build_relay_chart
%
%   Implements the sequence a real overcurrent + distance relay executes. This
%   is the dynamic counterpart of python/gridbench/protection.py: that module
%   answers "would the element operate for this fault", this chart answers
%   "what does the scheme actually DO over time", including the reclose cycle
%   and lockout that a single-shot calculation cannot express.
%
%   States
%     MONITOR    quiescent; watching current and apparent impedance
%     PICKUP     threshold exceeded, inverse-time integration running. A relay
%                does NOT trip on pickup -- it must stay picked up long enough,
%                which is what gives coordination its selectivity.
%     TRIP       trip signal asserted to the breaker
%     OPEN       breaker open, dead time running
%     RECLOSE    breaker reclosed; most faults are transient and clear here
%     LOCKOUT    permanent fault, reclose attempts exhausted, stays open
%
%   The reclose cycle is why this belongs in Stateflow: it is a timed, counted,
%   latching sequence, and the count of attempts before lockout is exactly the
%   kind of state a block diagram cannot hold cleanly.
%
%   THE IBR CONNECTION. The WP4 study found the dominant failure mode is not a
%   mis-timed trip but NO PICKUP AT ALL -- inverter-limited fault current never
%   reaches the threshold. In this chart that shows up as never leaving MONITOR.
%   Drive it with the fault currents from the IBR-aware solver and watch it sit
%   still while the fault persists.

if nargin < 1 || isempty(modelName), modelName = 'gb_relay'; end

cfg  = gb_config();
prot = cfg.protection;

if bdIsLoaded(modelName), close_system(modelName, 0); end
sfnew(modelName);
set_param(modelName, 'StopTime', num2str(cfg.dynamics.sim_time_s));

root  = sfroot;
chart = root.find('-isa', 'Stateflow.Chart', '-and', 'Path', modelName);
chart.Name = 'Relay';
chart.ActionLanguage = 'C';

%% ------------------------------------------------------------------ I/O
addData(chart, 'Imag',     'Input',  'double');   % measured current, pu
addData(chart, 'Zapp',     'Input',  'double');   % apparent impedance, pu
addData(chart, 'trip',     'Output', 'double');
addData(chart, 'breaker',  'Output', 'double');   % 1 = closed, 0 = open
addData(chart, 'state_id', 'Output', 'double');
addData(chart, 'zone',     'Output', 'double');
addLocal(chart, 'integ');                          % inverse-time accumulator
addLocal(chart, 'attempts');

addConst(chart, 'I_PICKUP',  1.25);               % pickup_multiplier x rating
addConst(chart, 'TMS',       0.10);
addConst(chart, 'K_IEC',     0.14);               % IEC standard inverse
addConst(chart, 'ALPHA_IEC', 0.02);
addConst(chart, 'Z1_REACH',  prot.distance_relay.zone1_reach_pct / 100);
addConst(chart, 'Z2_REACH',  prot.distance_relay.zone2_reach_pct / 100);
addConst(chart, 'T_Z2',      prot.distance_relay.zone2_delay_s);
addConst(chart, 'T_DEAD',    0.50);               % reclose dead time
addConst(chart, 'MAX_SHOTS', 2);

%% --------------------------------------------------------------- states
sMonitor = addState(chart, 'MONITOR', [40 40 240 110], sprintf([ ...
    'entry:\n' ...
    'trip = 0; breaker = 1; state_id = 0; zone = 0;\n' ...
    'integ = 0;']));

% Inverse-time integration: the element trips when the accumulated
% (M^alpha - 1)/K product reaches 1, which is the IEC characteristic in
% integral form -- the correct way to handle a current that varies during
% the fault, which a single t = f(I) lookup cannot.
sPickup = addState(chart, 'PICKUP', [340 40 300 170], sprintf([ ...
    'entry:\n' ...
    'state_id = 1;\n' ...
    'during:\n' ...
    'if (Zapp <= Z1_REACH) { zone = 1; }\n' ...
    'else if (Zapp <= Z2_REACH) { zone = 2; }\n' ...
    'else { zone = 3; }\n' ...
    'M = Imag / I_PICKUP;\n' ...
    'if (M > 1.0) {\n' ...
    '  integ = integ + 0.001 * (pow(M, ALPHA_IEC) - 1.0) / (TMS * K_IEC);\n' ...
    '}']));
addLocal(chart, 'M');

sTrip = addState(chart, 'TRIP', [700 40 220 110], sprintf([ ...
    'entry:\n' ...
    'trip = 1; state_id = 2;\n' ...
    'attempts = attempts + 1;']));

sOpen = addState(chart, 'OPEN', [700 220 220 110], sprintf([ ...
    'entry:\n' ...
    'breaker = 0; trip = 0; state_id = 3;']));

sReclose = addState(chart, 'RECLOSE', [340 260 240 110], sprintf([ ...
    'entry:\n' ...
    'breaker = 1; state_id = 4;\n' ...
    'integ = 0;']));

sLockout = addState(chart, 'LOCKOUT', [40 260 240 110], sprintf([ ...
    'entry:\n' ...
    'breaker = 0; trip = 0; state_id = 5;']));

addDefault(chart, sMonitor);

%% ---------------------------------------------------------- transitions
% Pickup and dropout. Dropout is what makes the element secure against
% transient inrush -- without it, any brief excursion would start a trip.
addTrans(chart, sMonitor, sPickup,  '[Imag > I_PICKUP]');
addTrans(chart, sPickup,  sMonitor, '[Imag <= I_PICKUP]');

% Zone 1 is instantaneous; zone 2 is time-delayed for coordination;
% the inverse-time element trips when its integral completes.
addTrans(chart, sPickup, sTrip, '[zone == 1]');
addTrans(chart, sPickup, sTrip, 'after(T_Z2, sec)[zone == 2]');
addTrans(chart, sPickup, sTrip, '[integ >= 1.0]');

addTrans(chart, sTrip,    sOpen,    'after(0.05, sec)');   % breaker opening time
addTrans(chart, sOpen,    sReclose, 'after(T_DEAD, sec)[attempts < MAX_SHOTS]');
addTrans(chart, sOpen,    sLockout, '[attempts >= MAX_SHOTS]');
addTrans(chart, sReclose, sMonitor, 'after(0.10, sec)[Imag <= I_PICKUP]');
addTrans(chart, sReclose, sPickup,  '[Imag > I_PICKUP]');

Simulink.BlockDiagram.arrangeSystem(modelName);
save_system(modelName, fullfile(gb_root(), 'matlab', 'stateflow', [modelName '.slx']));

fprintf('Built relay chart: %s\n', modelName);
fprintf('  pickup %.2f pu, Z1 reach %.0f%%, Z2 delay %.2f s, %d reclose shots\n', ...
    1.25, prot.distance_relay.zone1_reach_pct, ...
    prot.distance_relay.zone2_delay_s, 2);
fprintf('\nVerify three behaviours:\n');
fprintf('  1. Imag = 3 pu, Zapp = 0.5  -> zone 1, immediate trip, then reclose\n');
fprintf('  2. Imag = 3 pu sustained    -> two shots, then LOCKOUT (state_id 5)\n');
fprintf('  3. Imag = 1.1 pu            -> NEVER leaves MONITOR. This is the\n');
fprintf('     IBR failure mode from WP4: inverter-limited fault current below\n');
fprintf('     pickup means the relay simply does not see the fault.\n');
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
