function [Ybus, Yf, Yt] = gb_ybus(mpc)
%GB_YBUS  Bus admittance matrix built from first principles.
%
%   [Ybus, Yf, Yt] = GB_YBUS(mpc)
%
%   Deliberately NOT a call to MATPOWER's makeYbus. This is the MATLAB twin of
%   python/gridbench/ybus.py, and gb_validate_ybus compares it against makeYbus
%   so that the "from scratch" claim is checkable rather than asserted.
%
%   Branch model (MATPOWER convention), from bus f to bus t:
%       ys  = 1/(r + jx)                series admittance
%       t   = tau * exp(j*theta)        complex tap, tap on the FROM side
%       Yff = (ys + j*b/2) / tau^2
%       Yft = -ys / conj(t)
%       Ytf = -ys / t
%       Ytt =  ys + j*b/2

define_constants;

nb = size(mpc.bus, 1);
nl = size(mpc.branch, 1);
baseMVA = mpc.baseMVA;

% external bus numbers -> internal indices
e2i = sparse(mpc.bus(:, BUS_I), 1, 1:nb);
f = full(e2i(mpc.branch(:, F_BUS)));
t = full(e2i(mpc.branch(:, T_BUS)));

status = mpc.branch(:, BR_STATUS) > 0;

z  = mpc.branch(:, BR_R) + 1j * mpc.branch(:, BR_X);
ys = zeros(nl, 1);
ok = status & (abs(z) > 0);
ys(ok) = 1 ./ z(ok);

bc = mpc.branch(:, BR_B) .* status;

tau = mpc.branch(:, TAP);
tau(tau == 0) = 1;                       % TAP == 0 means nominal ratio
tp = tau .* exp(1j * pi/180 * mpc.branch(:, SHIFT));

Yff = (ys + 1j * bc / 2) ./ (tau .^ 2);
Yft = -ys ./ conj(tp);
Ytf = -ys ./ tp;
Ytt =  ys + 1j * bc / 2;

% branch-end admittance matrices (useful for flows)
i_idx = [(1:nl)'; (1:nl)'];
Yf = sparse(i_idx, [f; t], [Yff; Yft], nl, nb);
Yt = sparse(i_idx, [f; t], [Ytf; Ytt], nl, nb);

% assemble the bus matrix
Ybus = sparse([f; f; t; t], [f; t; f; t], [Yff; Yft; Ytf; Ytt], nb, nb);

% fixed bus shunts, specified as MW/MVAr demanded at V = 1.0 pu
Ysh = (mpc.bus(:, GS) + 1j * mpc.bus(:, BS)) / baseMVA;
Ybus = Ybus + sparse(1:nb, 1:nb, Ysh, nb, nb);
end
