function res = gb_fdlf(mpc, variant, opt)
%GB_FDLF  Fast-decoupled load flow, XB or BX formulation.
%
%   res = GB_FDLF(mpc)              % defaults to 'fdxb'
%   res = GB_FDLF(mpc, 'fdbx')
%
%   MATLAB twin of python/gridbench/solvers.py:fast_decoupled.
%
%   Exploits the weak P-|V| / Q-theta coupling of transmission networks (high
%   X/R). B' and B'' are constant, so they are factorised ONCE and each
%   iteration is far cheaper than a Newton step -- at the cost of linear rather
%   than quadratic convergence.
%
%   Worth stating in the report: the high-X/R assumption is what makes this
%   work, and it is exactly the assumption that fails on distribution feeders
%   and on converter-dominated networks with virtual impedance.
%
%   B' : line charging and taps removed; resistance also removed for XB.
%   B'': phase shifters removed;         resistance also removed for BX.

define_constants;
if nargin < 2 || isempty(variant), variant = 'fdxb'; end
if nargin < 3, opt = struct(); end
if ~isfield(opt, 'tol'),     opt.tol = 1e-10;  end
if ~isfield(opt, 'maxIter'), opt.maxIter = 100; end

tic;
nb = size(mpc.bus, 1);
baseMVA = mpc.baseMVA;
Ybus = gb_ybus(mpc);

e2i  = sparse(mpc.bus(:, BUS_I), 1, 1:nb);
gbus = full(e2i(mpc.gen(:, GEN_BUS)));
on   = mpc.gen(:, GEN_STATUS) > 0;

Sbus = -(mpc.bus(:, PD) + 1j * mpc.bus(:, QD)) / baseMVA;
for k = find(on)'
    Sbus(gbus(k)) = Sbus(gbus(k)) + (mpc.gen(k, PG) + 1j * mpc.gen(k, QG)) / baseMVA;
end

V = mpc.bus(:, VM) .* exp(1j * pi/180 * mpc.bus(:, VA));
for k = find(on)'
    V(gbus(k)) = mpc.gen(k, VG) * exp(1j * angle(V(gbus(k))));
end

pv   = find(mpc.bus(:, BUS_TYPE) == PV);
pq   = find(mpc.bus(:, BUS_TYPE) == PQ);
pvpq = sort([pv; pq]);

% --- B' -----------------------------------------------------------------
m1 = mpc;
m1.branch(:, BR_B) = 0;
m1.branch(:, TAP)  = 1;
m1.bus(:, BS)      = 0;
if strcmpi(variant, 'fdxb'), m1.branch(:, BR_R) = 0; end
Bp = -imag(gb_ybus(m1));

% --- B'' ----------------------------------------------------------------
m2 = mpc;
m2.branch(:, SHIFT) = 0;
if strcmpi(variant, 'fdbx'), m2.branch(:, BR_R) = 0; end
Bpp = -imag(gb_ybus(m2));

% Factorise once. decomposition() keeps the factors and applies the correct
% permutations internally, so each iteration is a cheap solve without any
% hand-rolled permutation bookkeeping to get wrong.
Bp_red  = decomposition(Bp(pvpq, pvpq),  'lu');
Bpp_red = decomposition(Bpp(pq, pq),     'lu');

history = [];
converged = false;
it = 0;

while it < opt.maxIter
    it = it + 1;

    % --- P-theta half ---
    mis = (V .* conj(Ybus * V) - Sbus) ./ abs(V);
    dP = real(mis(pvpq));
    if ~isempty(dP)
        dTh = Bp_red \ (-dP);
        Va = angle(V); Va(pvpq) = Va(pvpq) + dTh;
        V = abs(V) .* exp(1j * Va);
    end

    % --- Q-V half ---
    mis = (V .* conj(Ybus * V) - Sbus) ./ abs(V);
    dQ = imag(mis(pq));
    if ~isempty(dQ)
        dVm = Bpp_red \ (-dQ);
        Vm = abs(V); Vm(pq) = Vm(pq) + dVm;
        V = Vm .* exp(1j * angle(V));
    end

    mis = V .* conj(Ybus * V) - Sbus;
    F = [real(mis(pvpq)); imag(mis(pq))];
    history(end+1) = max(abs(F)); %#ok<AGROW>
    if isempty(F) || max(abs(F)) < opt.tol
        converged = true;
        break
    end
end

res = struct( ...
    'algorithm',  lower(variant), ...
    'converged',  converged, ...
    'iterations', it, ...
    'V',          V, ...
    'Vm',         abs(V), ...
    'Va_deg',     180/pi * angle(V), ...
    'mismatch',   history(end), ...
    'history',    history, ...
    'elapsed_s',  toc, ...
    'jacobian_cond', NaN, ...
    'q_limits_hit',  []);
end
