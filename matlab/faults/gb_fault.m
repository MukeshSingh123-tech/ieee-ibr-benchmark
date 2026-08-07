function res = gb_fault(seq, busNo, faultType, method, opt)
%GB_FAULT  Short-circuit solve: classical symmetrical components, or IBR-aware.
%
%   res = GB_FAULT(seq, busNo, '3LG')                  % classical
%   res = GB_FAULT(seq, busNo, 'SLG', 'ibr_aware')
%   res = GB_FAULT(seq, busNo, 'LL', 'ibr_aware', struct('k2',4,'iLimit',1.2))
%
%   MATLAB twin of python/gridbench/faults.py.
%
%   CLASSICAL: machines are voltage sources behind X"d, the three sequence
%   networks are independent, superposition holds, closed form.
%
%   IBR-AWARE: grid-following inverters are current-limited controlled sources
%   whose output depends on the very terminal voltage the fault produces, so
%   the network is NONLINEAR and must be iterated. Negative-sequence current
%   becomes a control choice (IEEE Std 2800-2022: I2 = K2*V2) rather than a
%   machine property, so the sequence networks stop being independent.
%
%   Both paths share one boundary-condition routine written in terms of
%   OPEN-CIRCUIT sequence voltages at the fault bus. The classical case is
%   simply the special case V2_oc = V0_oc = 0, which is what lets the two be
%   compared without a second implementation.

if nargin < 4 || isempty(method), method = 'classical'; end
if nargin < 5, opt = struct(); end
cfg = gb_config();
if ~isfield(opt, 'k2'),     opt.k2     = cfg.ibr.ieee2800.k2;      end
if ~isfield(opt, 'iLimit'), opt.iLimit = cfg.ibr.current_limit_pu; end
if ~isfield(opt, 'zf'),     opt.zf     = 0;    end
if ~isfield(opt, 'zg'),     opt.zg     = 0;    end
if ~isfield(opt, 'tol'),    opt.tol    = 1e-8; end
if ~isfield(opt, 'maxIter'),opt.maxIter= 200;  end
if ~isfield(opt, 'relax'),  opt.relax  = 0.7;  end

define_constants;
k  = full(seq.e2i(busNo));
n  = size(seq.Y1, 1);
z1 = seq.Z1(k, k); z2 = seq.Z2(k, k); z0 = seq.Z0(k, k);
ek = zeros(n, 1); ek(k) = 1;

%% ------------------------------------------------------------- classical
if strcmpi(method, 'classical')
    v1oc = seq.Vprefault(k);
    [i1, i2, i0] = faultCurrents(faultType, v1oc, 0, 0, z1, z2, z0, opt.zf, opt.zg);
    res = pack(faultType, busNo, 'classical', i1, i2, i0, ...
        seq.Vprefault - seq.Z1 * (ek * i1), ...
        -(seq.Z2 * (ek * i2)), -(seq.Z0 * (ek * i0)), 1, true, true, '');
    return
end

%% ------------------------------------------------------------- IBR-aware
% Well-posedness needs at least one VOLTAGE source. Synchronous machines
% provide one; so do grid-forming converters. With neither, every source is a
% current injection whose magnitude depends on a voltage nothing establishes.
if isempty(seq.voltageSourceBuses)
    res = pack(faultType, busNo, 'ibr_aware', 0, 0, 0, ...
        seq.Vprefault, zeros(n,1), zeros(n,1), 0, false, false, ...
        ['no voltage source remains (no synchronous plant, no grid-forming ' ...
         'converters) -- phasor fault model is ill-posed; add GFM or run EMT']);
    return
end

gfl = seq.gflBuses;
mpc = seq.mpc;
on  = mpc.gen(:, GEN_STATUS) > 0;
genBus = full(seq.e2i(mpc.gen(:, GEN_BUS)));

Psched = zeros(n, 1);
for i = gfl(:)'
    Psched(i) = sum(mpc.gen(on & genBus == i, PG)) / mpc.baseMVA;
end

V1sync = seq.Z1 * seq.e_internal;      % linear, fixed across iterations

inj1 = zeros(n, 1); inj2 = zeros(n, 1);
V1 = seq.Vprefault; V2 = zeros(n, 1); V0 = zeros(n, 1);
i1 = 0; i2 = 0; i0 = 0;
converged = false;
prevDelta = inf;
relax = opt.relax;
it = 0;

% With no grid-following units the problem is LINEAR -- the grid-forming and
% synchronous sources are already in the admittance matrix, so a single solve
% is exact and no iteration is required.
if isempty(gfl)
    v1oc = V1sync(k);
    [i1, i2, i0] = faultCurrents(faultType, v1oc, 0, 0, z1, z2, z0, opt.zf, opt.zg);
    res = pack(faultType, busNo, 'ibr_aware', i1, i2, i0, ...
        V1sync - seq.Z1 * (ek * i1), -(seq.Z2 * (ek * i2)), ...
        -(seq.Z0 * (ek * i0)), 1, true, true, ...
        'all voltage-source (GFM/synchronous): linear, solved directly');
    return
end

while it < opt.maxIter
    it = it + 1;

    new1 = zeros(n, 1); new2 = zeros(n, 1);
    for i = gfl(:)'
        [a, b] = ibrInjection(V1(i), V2(i), Psched(i), opt.iLimit, opt.k2, ...
                              seq.Vprefault(i));
        new1(i) = a; new2(i) = b;
    end
    inj1 = inj1 + relax * (new1 - inj1);
    inj2 = inj2 + relax * (new2 - inj2);

    V1ocVec = V1sync + seq.Z1 * inj1;
    V2ocVec = seq.Z2 * inj2;
    [i1, i2, i0] = faultCurrents(faultType, V1ocVec(k), V2ocVec(k), 0, ...
                                 z1, z2, z0, opt.zf, opt.zg);

    V1new = V1ocVec - seq.Z1 * (ek * i1);
    V2new = V2ocVec - seq.Z2 * (ek * i2);
    V0    = -(seq.Z0 * (ek * i0));

    delta = max([max(abs(V1new - V1)), max(abs(V2new - V2))]);
    V1 = V1new; V2 = V2new;

    % adaptive under-relaxation: the injection is a hard-clamped function of
    % the voltage it produces, so a fixed step can limit-cycle at the clamp
    if delta > prevDelta, relax = max(0.05, relax / 2); end
    prevDelta = delta;

    if delta < opt.tol, converged = true; break; end
end

res = pack(faultType, busNo, 'ibr_aware', i1, i2, i0, V1, V2, V0, it, ...
           converged, true, '');
end


%% ======================================================================
function [i1, i2, i0] = faultCurrents(faultType, v1oc, v2oc, v0oc, z1, z2, z0, zf, zg)
%FAULTCURRENTS  Sequence fault currents from open-circuit sequence voltages.
%   Generalised to nonzero v2oc/v0oc, which is what IBR control injection
%   produces. Setting v2oc = v0oc = 0 recovers every textbook formula exactly.
za = z1 + zf; zb = z2 + zf; zc = z0 + zf + 3 * zg;

switch upper(faultType)
    case '3LG'      % three phases tied together: V1' = V2' = 0, no ground path
        i1 = v1oc / za;  i2 = v2oc / zb;  i0 = 0;

    case 'SLG'      % phase a to ground: I1 = I2 = I0
        i = (v1oc + v2oc + v0oc) / (z1 + z2 + z0 + 3*zf + 3*zg);
        i1 = i; i2 = i; i0 = i;

    case 'LL'       % phases b-c, no ground: I0 = 0, I2 = -I1
        i1 = (v1oc - v2oc) / (z1 + z2 + zf);
        i2 = -i1; i0 = 0;

    case 'LLG'      % phases b-c to ground: I1+I2+I0 = 0, V1' = V2' = V0'
        A = [  1,   1,   1;
             -za,  zb,   0;
             -za,   0,  zc];
        rhs = [0; v2oc - v1oc; v0oc - v1oc];
        x = A \ rhs;
        i1 = x(1); i2 = x(2); i0 = x(3);

    otherwise
        error('gb_fault:type', 'unknown fault type %s', faultType);
end
end


function [i1, i2] = ibrInjection(v1, v2, pPU, iLimit, k2, vPre)
%IBRINJECTION  Current a grid-following inverter injects during a fault.
%   Reactive priority per IEEE 2800: Iq = K1*(1-|V|) takes the headroom first,
%   active power gets what is left. Negative sequence I2 = -K2*V2 is the clause
%   that keeps protection relays able to see the fault at all.
vm = abs(v1);
k1 = 2.0;
iq = min(k1 * max(0, 1 - vm), iLimit);
headroom = sqrt(max(0, iLimit^2 - iq^2));
if vm > 0.05
    ip = min(pPU / vm, headroom);
else
    ip = headroom;                       % saturated, continuous as vm -> 0
end

% PLL angle reference, with coast-through on voltage collapse. Below ~0.1 pu a
% real PLL loses lock and holds its last angle; tracking V/|V| there is
% numerically meaningless and destroys convergence.
if abs(vPre) > 1e-12, refHold = vPre / abs(vPre); else, refHold = 1; end
vLo = 0.02; vHi = 0.10;
if vm >= vHi
    ref = v1 / vm;
elseif vm <= vLo
    ref = refHold;
else
    w = (vm - vLo) / (vHi - vLo);
    blend = w * (v1 / vm) + (1 - w) * refHold;
    if abs(blend) > 1e-12, ref = blend / abs(blend); else, ref = refHold; end
end

i1 = ref * complex(ip, iq);
i2 = -k2 * v2;
if abs(i2) > iLimit, i2 = i2 / abs(i2) * iLimit; end
end


function res = pack(ft, bus, method, i1, i2, i0, V1, V2, V0, it, conv, wp, note)
a = exp(2j*pi/3);
A = [1 1 1; 1 a^2 a; 1 a a^2];
iph = A * [i0; i1; i2];
[~, worst] = max(abs(iph));

if abs(i1) > 1e-12, negRatio = abs(i2) / abs(i1); else, negRatio = NaN; end

res = struct('faultType', ft, 'bus', bus, 'method', method, ...
    'i1', i1, 'i2', i2, 'i0', i0, 'iPhase', iph, ...
    'iFaultMag', abs(iph(worst)), ...
    'iFaultAngleDeg', gb_wrap_deg(180/pi * angle(iph(worst))), ...
    'negSeqRatio', negRatio, ...
    'V1', V1, 'V2', V2, 'V0', V0, ...
    'iterations', it, 'converged', conv, 'wellPosed', wp, 'note', note);
end
