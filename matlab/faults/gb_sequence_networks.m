function seq = gb_sequence_networks(mpc, Vpre, ibrBuses, gfmBuses)
%GB_SEQUENCE_NETWORKS  Build positive/negative/zero sequence networks for a fault study.
%
%   seq = GB_SEQUENCE_NETWORKS(mpc, Vpre)
%   seq = GB_SEQUENCE_NETWORKS(mpc, Vpre, ibrBuses)
%   seq = GB_SEQUENCE_NETWORKS(mpc, Vpre, ibrBuses, gfmBuses)
%
%   MATLAB twin of python/gridbench/faults.py:build_sequence_networks.
%   `ibrBuses` and `gfmBuses` are EXTERNAL bus numbers; gfmBuses must be a
%   subset of ibrBuses.
%
%   Three source types, each treated as what it physically is:
%
%     synchronous machine   voltage source behind X"d        -> shunt + Norton source
%     GRID-FORMING (GFM)    voltage source behind Zvirtual   -> shunt + Norton source
%     GRID-FOLLOWING (GFL)  current source, no internal EMF  -> NO shunt
%
%   A grid-following inverter gets no shunt because it has no internal voltage
%   behind an impedance -- it synchronises to the grid through a PLL. Giving it
%   one would silently restore the very assumption this project sets out to test.
%
%   The presence of ANY voltage source is what makes the fault problem
%   well-posed. That is why grid-forming converters restore solvability to a
%   network with no synchronous plant left (see gb_fault_ibr).
%
%   IMPORTANT: the IEEE cases carry no machine or zero-sequence data. Everything
%   used here is declared in config/scenarios.yaml (machine_data, sequence_data)
%   rather than buried in this file.

define_constants;
if nargin < 3, ibrBuses = []; end
if nargin < 4, gfmBuses = []; end

cfg = gb_config();
md  = cfg.machine_data;
sd  = cfg.sequence_data;

mpc = ext2int(mpc);
nb  = size(mpc.bus, 1);
e2i = sparse(mpc.bus(:, BUS_I), 1, 1:nb);

ibrIdx = internalIdx(e2i, ibrBuses);
gfmIdx = intersect(internalIdx(e2i, gfmBuses), ibrIdx);
gflIdx = setdiff(ibrIdx, gfmIdx);

on       = mpc.gen(:, GEN_STATUS) > 0;
genBus   = full(e2i(mpc.gen(:, GEN_BUS)));
allGen   = unique(genBus(on));
syncIdx  = setdiff(allGen, ibrIdx);

%% --- machine shunts ------------------------------------------------------
Ysh1 = machineShunt(mpc, md.xdpp_pu, md.r_over_x, syncIdx, genBus, on, cfg);
Ysh2 = machineShunt(mpc, md.x2_pu,   md.r_over_x, syncIdx, genBus, on, cfg);
x0eff = md.x0_pu + 3 * md.xn_pu;      % neutral impedance appears as 3*Zn
Ysh0 = machineShunt(mpc, x0eff, md.r_over_x, syncIdx, genBus, on, cfg);

%% --- grid-forming converters: voltage source behind virtual impedance ----
Ygfm = zeros(nb, 1);
zv = cfg.dynamics.gfm.virtual_impedance_pu;
for i = gfmIdx(:)'
    rows = find(on & genBus == i);
    if isempty(rows), continue; end
    Srat = gb_inverter_rating(mpc, rows);          % pu on system base
    if Srat <= 0, continue; end
    zSys = complex(zv.r, zv.x) / Srat;             % converter base -> system base
    Ygfm(i) = Ygfm(i) + 1 / zSys;
end

%% --- loads as constant impedance from the prefault solution --------------
Yload = zeros(nb, 1);
if ~strcmpi(sd.load_model, 'neglect')
    Sload = (mpc.bus(:, PD) + 1j * mpc.bus(:, QD)) / mpc.baseMVA;
    vm2 = abs(Vpre) .^ 2; vm2(vm2 < 1e-9) = 1;
    Yload = conj(Sload) ./ vm2;
end

%% --- assemble -----------------------------------------------------------
Y1 = gb_ybus(mpc) + sparse(1:nb, 1:nb, Ysh1 + Ygfm       + Yload, nb, nb);
% GFM control suppresses negative-sequence current -> higher effective Z (2x)
Y2 = gb_ybus(mpc) + sparse(1:nb, 1:nb, Ysh2 + Ygfm / 2   + Yload, nb, nb);

% zero sequence: synthesised from positive-sequence data
m0 = mpc;
m0.branch(:, BR_R) = m0.branch(:, BR_R) * sd.line_r0_over_r1;
m0.branch(:, BR_X) = m0.branch(:, BR_X) * sd.line_x0_over_x1;
m0.branch(:, BR_B) = m0.branch(:, BR_B) * sd.line_b0_over_b1;
if sd.transformer_blocks_zero_seq
    % delta-wye(g) blocks zero sequence; transformers are branches with a tap
    isXfmr = m0.branch(:, TAP) ~= 0 & m0.branch(:, TAP) ~= 1;
    m0.branch(isXfmr, BR_STATUS) = 0;
end
% converters are interfaced through delta-wye, so contribute no zero sequence
Y0 = gb_ybus(m0) + sparse(1:nb, 1:nb, Ysh0 + 1e-9, nb, nb);

%% --- Norton sources for every VOLTAGE source ----------------------------
Eint = zeros(nb, 1);
Eint(syncIdx) = Vpre(syncIdx) .* Ysh1(syncIdx);
if ~isempty(gfmIdx)
    Eint(gfmIdx) = Vpre(gfmIdx) .* Ygfm(gfmIdx);
end

seq = struct( ...
    'Y1', Y1, 'Y2', Y2, 'Y0', Y0, ...
    'Z1', inv(full(Y1)), 'Z2', inv(full(Y2)), 'Z0', inv(full(Y0)), ...
    'e_internal', Eint, ...
    'syncBuses', syncIdx(:), 'ibrBuses', ibrIdx(:), ...
    'gfmBuses', gfmIdx(:),  'gflBuses', gflIdx(:), ...
    'voltageSourceBuses', union(syncIdx(:), gfmIdx(:)), ...
    'Vprefault', Vpre, 'mpc', mpc, 'e2i', e2i);
end


%% ======================================================================
function idx = internalIdx(e2i, busNumbers)
if isempty(busNumbers)
    idx = [];
else
    idx = full(e2i(busNumbers(:)));
    idx = idx(idx > 0);
end
end


function Ysh = machineShunt(mpc, xpu, rOverX, buses, genBus, on, cfg)
%MACHINESHUNT  Diagonal shunt from machine reactance, referred to SYSTEM base.
%   The slack is an unbounded external grid with no meaningful machine rating,
%   so it gets a Thevenin impedance derived from an assumed short-circuit level
%   instead -- which is how a utility actually specifies an interconnection.
define_constants;
nb  = size(mpc.bus, 1);
Ysh = zeros(nb, 1);
refBuses = find(mpc.bus(:, BUS_TYPE) == REF);

scMVA = cfg.machine_data.slack.sc_mva;
if isempty(scMVA) || (isnumeric(scMVA) && isscalar(scMVA) && scMVA == 0)
    scMVA = cfg.machine_data.slack.sc_multiple_of_load * sum(mpc.bus(:, PD));
end
scPU = scMVA / mpc.baseMVA;
xr   = cfg.machine_data.slack.x_over_r;
zmag = 1 / scPU;
xs   = zmag * xr / sqrt(1 + xr^2);
Zslack = complex(xs / xr, xs);

for i = buses(:)'
    rows = find(on & genBus == i);
    for k = rows(:)'
        if ismember(i, refBuses)
            Ysh(i) = Ysh(i) + 1 / Zslack;
            continue
        end
        mbase = mpc.gen(k, PMAX);
        if ~isfinite(mbase) || mbase <= 0 || mbase >= 1e6
            mbase = abs(mpc.gen(k, PG));
        end
        if mbase <= 0, mbase = mpc.baseMVA; end
        xSys = xpu * mpc.baseMVA / mbase;
        Ysh(i) = Ysh(i) + 1 / complex(rOverX * xSys, xSys);
    end
end
end
