function Yred = gb_reduce(net, faultBus, outagedBranch)
%GB_REDUCE  Kron-reduce the network to the internal nodes of the dynamic units.
%
%   Yred = GB_REDUCE(net)                       % pre-fault
%   Yred = GB_REDUCE(net, faultBusIdx)          % during a bolted fault
%   Yred = GB_REDUCE(net, [], branchIdx)        % post-fault, branch tripped
%
%   Loads, and grid-following converters, become fixed shunt admittances from
%   the prefault solution. A grid-following unit has no rotor and no internal
%   EMF, so a constant-impedance representation is the consistent choice: it can
%   neither absorb nor release rotational energy.

define_constants;
if nargin < 2, faultBus = []; end
if nargin < 3, outagedBranch = []; end

mpc = net.mpc;
if ~isempty(outagedBranch)
    mpc.branch(outagedBranch, BR_STATUS) = 0;
end

nb = size(mpc.bus, 1);
Y = full(gb_ybus(mpc));

% loads as constant impedance
Sload = (mpc.bus(:, PD) + 1j * mpc.bus(:, QD)) / mpc.baseMVA;
vm2 = abs(net.Vpre) .^ 2; vm2(vm2 < 1e-9) = 1;
Y = Y + diag(conj(Sload) ./ vm2);

% grid-following injections as constant impedance (negative load)
on = mpc.gen(:, GEN_STATUS) > 0;
genBus = full(net.e2i(mpc.gen(:, GEN_BUS)));
for i = net.gflBuses(:)'
    rows = find(on & genBus == i);
    if isempty(rows), continue; end
    P = sum(mpc.gen(rows, PG)) / mpc.baseMVA;
    Y(i, i) = Y(i, i) + conj(-complex(P, 0)) / vm2(i);
end

% bolted three-phase fault
if ~isempty(faultBus)
    Y(faultBus, faultBus) = Y(faultBus, faultBus) + 1e7;
end

% augment with the internal nodes behind X'd
nDyn = net.n;
Yaug = zeros(nb + nDyn);
Yaug(1:nb, 1:nb) = Y;
for k = 1:nDyn
    i = net.dynBuses(k); node = nb + k;
    ys = 1 / (1j * net.Xtr(k));
    Yaug(node, node) = Yaug(node, node) + ys;
    Yaug(i, i)       = Yaug(i, i)       + ys;
    Yaug(node, i)    = Yaug(node, i)    - ys;
    Yaug(i, node)    = Yaug(i, node)    - ys;
end

keep = (nb+1):(nb+nDyn);
drop = 1:nb;
Ykk = Yaug(keep, keep);
Ykd = Yaug(keep, drop);
Ydk = Yaug(drop, keep);
Ydd = Yaug(drop, drop);
Yred = Ykk - Ykd * (Ydd \ Ydk);
end
