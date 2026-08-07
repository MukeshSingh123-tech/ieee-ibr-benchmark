function res = gb_cct(mpc, Vpre, net, faultBus, ibrGens, gfmBuses, outagedBranch, tMax)
%GB_CCT  Critical clearing time by bisection on the swing equations.
%
%   res = GB_CCT(mpc, Vpre, net, faultBus, ibrGens, gfmBuses, outagedBranch, tMax)
%
%   MATLAB twin of python/gridbench/transient.py:critical_clearing_time.
%
%   Integrates the classical multi-machine swing model through the fault-on and
%   post-fault periods with RK4, and bisects on clearing time. Instability is
%   declared when any rotor angle departs from the centre of inertia by more
%   than 180 degrees -- the standard loss-of-synchronism criterion.
%
%   See the header of wp7_transient_matlab for why the resulting CCT must not be
%   compared across scenarios with different numbers of dynamic units.

if nargin < 8 || isempty(tMax), tMax = 1.0; end
if nargin < 7, outagedBranch = []; end

fb = full(net.e2i(faultBus));
Yf = gb_reduce(net, fb, []);
Yp = gb_reduce(net, [], outagedBranch);

tol = 0.005;

if ~swingStable(net, Yf, Yp, 0.0)
    res = struct('cct', 0.0, 'bracketed', false, ...
                 'note', 'unstable even at zero clearing time');
    return
end
if swingStable(net, Yf, Yp, tMax)
    res = struct('cct', tMax, 'bracketed', false, ...
                 'note', sprintf('still stable at tMax=%.2fs', tMax));
    return
end

lo = 0.0; hi = tMax;
while (hi - lo) > tol
    mid = 0.5 * (lo + hi);
    if swingStable(net, Yf, Yp, mid), lo = mid; else, hi = mid; end
end
res = struct('cct', lo, 'bracketed', true, 'note', '');
end


function stable = swingStable(net, Yf, Yp, tClear)
%SWINGSTABLE  RK4 integration of the swing equations; loss-of-synchronism test.
f0 = 60; ws = 2*pi*f0;
dt = 0.002; tEnd = 5.0;
n = net.n;

delta = net.delta0; omega = zeros(n, 1);
Htot = sum(net.H);
stable = true;

for k = 1:round(tEnd/dt)
    t = (k-1)*dt;
    if t < tClear, Y = Yf; else, Y = Yp; end

    [k1d, k1w] = deriv(net, delta,               omega,               Y, ws);
    [k2d, k2w] = deriv(net, delta + 0.5*dt*k1d,  omega + 0.5*dt*k1w,  Y, ws);
    [k3d, k3w] = deriv(net, delta + 0.5*dt*k2d,  omega + 0.5*dt*k2w,  Y, ws);
    [k4d, k4w] = deriv(net, delta + dt*k3d,      omega + dt*k3w,      Y, ws);

    delta = delta + dt/6 * (k1d + 2*k2d + 2*k3d + k4d);
    omega = omega + dt/6 * (k1w + 2*k2w + 2*k3w + k4w);

    if ~all(isfinite(delta)), stable = false; return; end
    coi = sum(net.H .* delta) / Htot;
    if max(abs(delta - coi)) > pi
        stable = false; return
    end
end
end


function [dd, dw] = deriv(net, delta, omega, Y, ws)
V  = abs(net.E) .* exp(1j*delta);
Pe = real(V .* conj(Y * V));
dd = omega * ws;
dw = (net.Pm - Pe - net.D .* omega) ./ (2 * net.H);
end
