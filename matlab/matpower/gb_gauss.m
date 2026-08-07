function res = gb_gauss(mpc, opt)
%GB_GAUSS  Gauss-Seidel power flow with an acceleration factor.
%
%   res = GB_GAUSS(mpc)
%   res = GB_GAUSS(mpc, struct('tol',1e-8,'maxIter',8000,'accel',1.6))
%
%   MATLAB twin of python/gridbench/solvers.py:gauss_seidel.
%
%   For each PQ bus i:      V_i <- (conj(S_i)/conj(V_i) - sum_{k~=i} Y_ik V_k) / Y_ii
%   For each PV bus i:      recompute Q from the present voltage, respect the
%                           Q limits, then keep only the NEW ANGLE and restore
%                           the scheduled magnitude.
%
%   Convergence is linear and the iteration count scales badly with system size.
%   That degradation is the WP1 result, so maxIter is set high enough to let the
%   method finish rather than cutting it off early and calling it divergence.
%   The acceleration factor matters: IEEE 39-bus converges at accel <= 1.4 and
%   diverges at 1.6, which Newton-Raphson has no equivalent of.

define_constants;
if nargin < 2, opt = struct(); end
if ~isfield(opt, 'tol'),     opt.tol = 1e-8;     end
if ~isfield(opt, 'maxIter'), opt.maxIter = 8000; end
if ~isfield(opt, 'accel'),   opt.accel = 1.6;    end

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

ref = find(mpc.bus(:, BUS_TYPE) == REF);
pv  = find(mpc.bus(:, BUS_TYPE) == PV);

Vg = abs(V);
Qmax = zeros(nb, 1); Qmin = zeros(nb, 1);
for b = 1:nb
    rows = find(on & gbus == b);
    if ~isempty(rows)
        Qmax(b) = sum(mpc.gen(rows, QMAX)) / baseMVA;
        Qmin(b) = sum(mpc.gen(rows, QMIN)) / baseMVA;
    end
end

isRef = false(nb, 1); isRef(ref) = true;
isPV  = false(nb, 1); isPV(pv)   = true;

history = [];
converged = false;
it = 0;

while it < opt.maxIter
    it = it + 1;
    Vprev = V;

    for i = 1:nb
        if isRef(i), continue; end

        rowSum = Ybus(i, :) * V - Ybus(i, i) * V(i);

        if isPV(i)
            qCalc = -imag(conj(V(i)) * (rowSum + Ybus(i, i) * V(i)));
            qLoad = mpc.bus(i, QD) / baseMVA;
            qGen  = qCalc + qLoad;
            if qGen > Qmax(i)
                Si = complex(real(Sbus(i)), Qmax(i) - qLoad);
            elseif qGen < Qmin(i)
                Si = complex(real(Sbus(i)), Qmin(i) - qLoad);
            else
                Si = complex(real(Sbus(i)), qCalc);
                Vnew = (conj(Si) / conj(V(i)) - rowSum) / Ybus(i, i);
                V(i) = Vg(i) * exp(1j * angle(Vnew));   % hold scheduled |V|
                continue
            end
        else
            Si = Sbus(i);
        end

        Vnew = (conj(Si) / conj(V(i)) - rowSum) / Ybus(i, i);
        V(i) = V(i) + opt.accel * (Vnew - V(i));
    end

    delta = max(abs(V - Vprev));
    history(end+1) = delta; %#ok<AGROW>
    if delta < opt.tol
        converged = true;
        break
    end
end

res = struct( ...
    'algorithm',  'gs', ...
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
