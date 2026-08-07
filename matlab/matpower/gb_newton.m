function res = gb_newton(mpc, opt)
%GB_NEWTON  Newton-Raphson power flow written from first principles.
%
%   res = GB_NEWTON(mpc)
%   res = GB_NEWTON(mpc, struct('tol',1e-10,'maxIter',50,'enforceQ',true))
%
%   MATLAB twin of python/gridbench/solvers.py:newton_raphson. Returns the same
%   fields as the Python PFResult so the two can be compared column-by-column
%   by compare/gb_compare_python.m.
%
%   Polar form. Unknowns are theta at PV+PQ buses and |V| at PQ buses; the
%   mismatch is [dP(pv,pq); dQ(pq)] and the Jacobian is the usual 4-block
%
%       J = [ dP/dtheta   dP/d|V| ;
%             dQ/dtheta   dQ/d|V| ]
%
%   Convergence is quadratic: on a well-conditioned IEEE case this reaches
%   1e-10 in 3-5 iterations regardless of how many buses there are.

define_constants;
if nargin < 2, opt = struct(); end
if ~isfield(opt, 'tol'),      opt.tol = 1e-10;  end
if ~isfield(opt, 'maxIter'),  opt.maxIter = 50; end
if ~isfield(opt, 'enforceQ'), opt.enforceQ = false; end

tic;
nb = size(mpc.bus, 1);
baseMVA = mpc.baseMVA;
Ybus = gb_ybus(mpc);

e2i = sparse(mpc.bus(:, BUS_I), 1, 1:nb);
gbus = full(e2i(mpc.gen(:, GEN_BUS)));
on = mpc.gen(:, GEN_STATUS) > 0;

% scheduled injections
Sbus = -(mpc.bus(:, PD) + 1j * mpc.bus(:, QD)) / baseMVA;
for k = find(on)'
    Sbus(gbus(k)) = Sbus(gbus(k)) + (mpc.gen(k, PG) + 1j * mpc.gen(k, QG)) / baseMVA;
end

% flat start with generator voltage setpoints applied
V = mpc.bus(:, VM) .* exp(1j * pi/180 * mpc.bus(:, VA));
for k = find(on)'
    V(gbus(k)) = mpc.gen(k, VG) * exp(1j * angle(V(gbus(k))));
end

ref = find(mpc.bus(:, BUS_TYPE) == REF);
pv  = find(mpc.bus(:, BUS_TYPE) == PV);
pq  = find(mpc.bus(:, BUS_TYPE) == PQ);

history = [];
totalIter = 0;
switched = [];
converged = false;

for outer = 1:(1 + 9 * opt.enforceQ)
    pvpq = sort([pv; pq]);

    mis = V .* conj(Ybus * V) - Sbus;
    F = [real(mis(pvpq)); imag(mis(pq))];
    history(end+1) = max(abs(F)); %#ok<AGROW>
    converged = isempty(F) || max(abs(F)) < opt.tol;

    it = 0;
    while ~converged && it < opt.maxIter
        it = it + 1; totalIter = totalIter + 1;

        J = gb_jacobian(Ybus, V, pvpq, pq);
        dx = -(J \ F);

        Va = angle(V); Vm = abs(V);
        Va(pvpq) = Va(pvpq) + dx(1:numel(pvpq));
        Vm(pq)   = Vm(pq)   + dx(numel(pvpq)+1:end);
        V = Vm .* exp(1j * Va);

        mis = V .* conj(Ybus * V) - Sbus;
        F = [real(mis(pvpq)); imag(mis(pq))];
        history(end+1) = max(abs(F)); %#ok<AGROW>
        converged = max(abs(F)) < opt.tol;
    end

    if ~opt.enforceQ || ~converged, break; end

    % --- Q-limit enforcement: convert violating PV buses to PQ -----------
    violated = false;
    S = V .* conj(Ybus * V);
    for b = pv'
        rows = find(on & gbus == b);
        if isempty(rows), continue; end
        qgen = (imag(S(b)) + mpc.bus(b, QD) / baseMVA) * baseMVA;
        qmax = sum(mpc.gen(rows, QMAX));
        qmin = sum(mpc.gen(rows, QMIN));
        if qgen > qmax + 1e-6 || qgen < qmin - 1e-6
            qfix = min(max(qgen, qmin), qmax);
            Sbus(b) = Sbus(b) + 1j * (qfix - qgen) / baseMVA;
            pv = pv(pv ~= b);
            pq = sort([pq; b]);
            switched(end+1) = mpc.bus(b, BUS_I); %#ok<AGROW>
            violated = true;
        end
    end
    if ~violated, break; end
end

condJ = NaN;
if converged
    try
        condJ = cond(full(gb_jacobian(Ybus, V, sort([pv; pq]), pq)));
    catch
    end
end

res = struct( ...
    'algorithm',  'nr', ...
    'converged',  converged, ...
    'iterations', totalIter, ...
    'V',          V, ...
    'Vm',         abs(V), ...
    'Va_deg',     180/pi * angle(V), ...
    'mismatch',   history(end), ...
    'history',    history, ...
    'elapsed_s',  toc, ...
    'jacobian_cond', condJ, ...
    'q_limits_hit',  switched);
end


function J = gb_jacobian(Ybus, V, pvpq, pq)
%GB_JACOBIAN  Full polar Newton-Raphson Jacobian.
%   dS/d|V|   = diag(V./|V|) conj(diag(Ibus)) + diag(V) conj(Ybus diag(V./|V|))
%   dS/dtheta = j diag(V) conj(diag(Ibus) - Ybus diag(V))
n     = numel(V);
Ibus  = Ybus * V;
Vnorm = V ./ abs(V);
diagV     = sparse(1:n, 1:n, V,     n, n);
diagIbus  = sparse(1:n, 1:n, Ibus,  n, n);
diagVnorm = sparse(1:n, 1:n, Vnorm, n, n);

dS_dVm = diagVnorm * conj(diagIbus) + diagV * conj(Ybus * diagVnorm);
dS_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV);

J = [ real(dS_dVa(pvpq, pvpq))  real(dS_dVm(pvpq, pq));
      imag(dS_dVa(pq,   pvpq))  imag(dS_dVm(pq,   pq)) ];
end
