function net = gb_dynamic_network(mpc, Vpre, ibrGens, gfmBuses, hVirtual)
%GB_DYNAMIC_NETWORK  Classical transient-stability model, Kron-reduced.
%
%   net = GB_DYNAMIC_NETWORK(mpc, Vpre)                    % all synchronous
%   net = GB_DYNAMIC_NETWORK(mpc, Vpre, ibrGens, gfmBuses)
%   net = GB_DYNAMIC_NETWORK(mpc, Vpre, ibrGens, gfmBuses, hVirtual)
%
%   MATLAB twin of python/gridbench/transient.py:build_dynamic_network.
%
%   Dynamic (voltage-source) units are the synchronous machines PLUS the
%   grid-forming converters -- a GFM unit under virtual-synchronous-machine
%   control obeys the same swing equation, with a synthesised inertia constant.
%   Grid-following converters are EXCLUDED: no rotor, no angle to swing. They
%   fold into the network as constant-impedance injections.
%
%   Machine data comes from config/scenarios.yaml (machine_dynamics), which
%   carries the published values for the WSCC 9-bus and IEEE 39-bus benchmarks,
%   already referred to the 100 MVA system base. Using a generic inertia
%   constant instead gives a CCT several times too long.

define_constants;
if nargin < 3, ibrGens = []; end
if nargin < 4, gfmBuses = []; end

cfg = gb_config();
if nargin < 5 || isempty(hVirtual)
    hVirtual = cfg.dynamics.gfm.inertia_const_s;
end
hDefault = cfg.dynamics.h_default_s;
xdpDef   = cfg.machine_data.xdp_pu;
zv       = cfg.dynamics.gfm.virtual_impedance_pu;
baseMVA  = mpc.baseMVA;

nb  = size(mpc.bus, 1);
e2i = sparse(mpc.bus(:, BUS_I), 1, 1:nb);
on  = mpc.gen(:, GEN_STATUS) > 0;
genBus = full(e2i(mpc.gen(:, GEN_BUS)));

gfmIdx = [];
if ~isempty(gfmBuses), gfmIdx = full(e2i(gfmBuses(:))); end
ibrIdx = [];
if ~isempty(ibrGens),  ibrIdx = unique(genBus(ibrGens)); end
gflIdx = setdiff(ibrIdx, gfmIdx);

% published per-machine dynamic data, if this case is a benchmark
dyn = struct();
caseName = guessCaseName(mpc);
if ~isempty(caseName) && isfield(cfg, 'machine_dynamics') && ...
        isfield(cfg.machine_dynamics, caseName)
    dyn = cfg.machine_dynamics.(caseName);
end

%% --- which generator rows are dynamic ----------------------------------
gflGenRows = [];
if ~isempty(ibrGens)
    gflGenRows = ibrGens(~ismember(genBus(ibrGens), gfmIdx));
end
dynRows = find(on);
dynRows = dynRows(~ismember(dynRows, gflGenRows));

nDyn = numel(dynRows);
dynBuses = zeros(nDyn, 1); Xtr = zeros(nDyn, 1);
H = zeros(nDyn, 1); isGfm = false(nDyn, 1);

for k = 1:nDyn
    g = dynRows(k);
    b = mpc.gen(g, GEN_BUS);
    dynBuses(k) = full(e2i(b));

    rating = mpc.gen(g, PMAX);
    if ~isfinite(rating) || rating <= 0 || rating >= 1e6
        rating = abs(mpc.gen(g, PG));
    end
    if rating <= 0, rating = baseMVA; end

    fld = sprintf('x%d', b);          % jsondecode turns numeric keys into xN
    hasPub = isstruct(dyn) && isfield(dyn, fld);

    if ismember(full(e2i(b)), gfmIdx)
        Srat = gb_inverter_rating(mpc, g);
        if Srat > 0, Xtr(k) = zv.x / Srat; else, Xtr(k) = xdpDef; end
        H(k) = hVirtual * rating / baseMVA;
        isGfm(k) = true;
    elseif hasPub
        Xtr(k) = dyn.(fld).xdp;       % already on the system base
        H(k)   = dyn.(fld).h;
    else
        Xtr(k) = xdpDef * baseMVA / rating;
        H(k)   = hDefault * rating / baseMVA;
    end
end

%% --- internal EMFs from the prefault solution --------------------------
Ybus = gb_ybus(mpc);
S = Vpre .* conj(Ybus * Vpre);
E = zeros(nDyn, 1); Pm = zeros(nDyn, 1);
for k = 1:nDyn
    g = dynRows(k); i = dynBuses(k);
    Pg = mpc.gen(g, PG) / baseMVA;
    Qg = imag(S(i)) + mpc.bus(i, QD) / baseMVA;
    if abs(Vpre(i)) > 1e-9
        Ig = conj(complex(Pg, Qg) / Vpre(i));
    else
        Ig = 0;
    end
    E(k)  = Vpre(i) + 1j * Xtr(k) * Ig;
    Pm(k) = Pg;
end

net = struct('mpc', mpc, 'e2i', e2i, 'Vpre', Vpre, ...
    'dynRows', dynRows, 'dynBuses', dynBuses, 'Xtr', Xtr, ...
    'H', H, 'Pm', Pm, 'D', zeros(nDyn,1), 'E', E, ...
    'isGfm', isGfm, 'gflBuses', gflIdx, 'delta0', angle(E), 'n', nDyn);

net.Yred = gb_reduce(net, [], []);
end


function name = guessCaseName(mpc)
%GUESSCASENAME  Identify a standard benchmark by its size.
switch size(mpc.bus, 1)
    case 9,  name = 'case9';
    case 39, name = 'case39';
    otherwise, name = '';
end
end
