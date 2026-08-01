"""Weighted least-squares state estimation, bad-data detection, and false data
injection attacks.

This is the bridge from the power-systems work to the security work: the load
flow solutions from WP1 become the measurement set a control centre would
actually see, and the question becomes whether an attacker can corrupt the
operator's picture of the grid without being caught.

The classical defence is a chi-squared test on the measurement residual. Its
weakness is precise and well known (Liu, Ning & Reiter, 2009): if an attacker
knows the topology H, any attack vector of the form

        a = H c

adds exactly zero to the residual, because the estimator attributes it to a
genuine state change of c. Such an attack is UNDETECTABLE by residual testing at
any threshold -- the defence does not fail by being badly tuned, it fails
structurally.

That is the motivation for out-of-band verification (PMU physics checks,
cryptographic attestation of meter data), which is what the ChainPMU line of
work addresses. This module quantifies the gap that motivates it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np
from scipy import stats

from . import ybus as yb
from .ppc import BR_STATUS, F_BUS, T_BUS, bus_index, bus_types


MeasurementType = Literal["vm", "pinj", "qinj", "pflow", "qflow"]


# =============================================================================
# measurement model
# =============================================================================
@dataclass
class MeasurementSet:
    """A control-centre measurement set built from a solved power flow."""

    types: list[MeasurementType]
    indices: list[tuple[int, ...]]      # bus index, or (branch, from/to)
    values: np.ndarray                  # measured values (with noise)
    truth: np.ndarray                   # noise-free values, for scoring
    sigma: np.ndarray                   # standard deviations

    @property
    def n(self) -> int:
        return self.values.size

    @property
    def weights(self) -> np.ndarray:
        """W = diag(1/sigma^2) -- the inverse-covariance weighting of WLS."""
        return 1.0 / self.sigma ** 2


def _flows(ppc: dict[str, Any], v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sf, st = yb.branch_flows(ppc, v)
    return sf, st


def build_measurements(
    ppc: dict[str, Any],
    v: np.ndarray,
    sigma_v: float = 0.004,
    sigma_p: float = 0.008,
    redundancy: float = 1.0,
    seed: int = 20260801,
) -> MeasurementSet:
    """Construct a realistic, redundant measurement set from a solved power flow.

    Voltage magnitudes at every bus, injections at every bus, and flows at both
    ends of every in-service branch. Redundancy (measurements / states) of 2-4
    is typical in real EMS installations and is what makes bad-data detection
    possible at all -- with no redundancy there is nothing to cross-check.
    """
    rng = np.random.default_rng(seed)
    ybus = yb.build(ppc)
    s_bus = v * np.conj(ybus @ v)
    sf, st = _flows(ppc, v)

    types: list[MeasurementType] = []
    indices: list[tuple[int, ...]] = []
    truth: list[float] = []
    sigma: list[float] = []

    n = ppc["bus"].shape[0]
    for i in range(n):
        types.append("vm"); indices.append((i,))
        truth.append(float(abs(v[i]))); sigma.append(sigma_v)
    for i in range(n):
        types.append("pinj"); indices.append((i,))
        truth.append(float(s_bus[i].real)); sigma.append(sigma_p)
        types.append("qinj"); indices.append((i,))
        truth.append(float(s_bus[i].imag)); sigma.append(sigma_p)

    for k in range(ppc["branch"].shape[0]):
        if ppc["branch"][k, BR_STATUS] <= 0:
            continue
        if rng.random() > redundancy:
            continue
        types.append("pflow"); indices.append((k, 0))
        truth.append(float(sf[k].real)); sigma.append(sigma_p)
        types.append("qflow"); indices.append((k, 0))
        truth.append(float(sf[k].imag)); sigma.append(sigma_p)

    truth_arr = np.array(truth)
    sigma_arr = np.array(sigma)
    values = truth_arr + rng.normal(0.0, sigma_arr)
    return MeasurementSet(types, indices, values, truth_arr, sigma_arr)


def measurement_function(
    ppc: dict[str, Any], meas: MeasurementSet, v: np.ndarray,
) -> np.ndarray:
    """h(x) -- what the measurements WOULD read for a given state."""
    ybus = yb.build(ppc)
    s_bus = v * np.conj(ybus @ v)
    sf, st = _flows(ppc, v)

    out = np.empty(meas.n)
    for m, (kind, idx) in enumerate(zip(meas.types, meas.indices)):
        if kind == "vm":
            out[m] = abs(v[idx[0]])
        elif kind == "pinj":
            out[m] = s_bus[idx[0]].real
        elif kind == "qinj":
            out[m] = s_bus[idx[0]].imag
        elif kind == "pflow":
            out[m] = sf[idx[0]].real
        else:
            out[m] = sf[idx[0]].imag
    return out


def jacobian(
    ppc: dict[str, Any], meas: MeasurementSet, v: np.ndarray,
) -> np.ndarray:
    """Numerical measurement Jacobian H = dh/dx, x = [theta(non-ref); |V|].

    Numerical rather than analytic: the analytic forms for flow measurements
    through off-nominal taps and phase shifters are error-prone, and at IEEE
    test-system scale the cost is irrelevant. `test_estimation.py` checks the
    estimator recovers the true state, which would fail if H were wrong.
    """
    ref, pv, pq = bus_types(ppc)
    n = v.size
    non_ref = np.array([i for i in range(n) if i not in set(ref.tolist())])
    n_state = non_ref.size + n

    h0 = measurement_function(ppc, meas, v)
    jac = np.empty((meas.n, n_state))
    eps = 1e-7

    for j, i in enumerate(non_ref):
        vp = v.copy()
        vp[i] = abs(vp[i]) * np.exp(1j * (np.angle(vp[i]) + eps))
        jac[:, j] = (measurement_function(ppc, meas, vp) - h0) / eps

    for j in range(n):
        vp = v.copy()
        vp[j] = (abs(vp[j]) + eps) * np.exp(1j * np.angle(vp[j]))
        jac[:, non_ref.size + j] = (measurement_function(ppc, meas, vp) - h0) / eps

    return jac


# =============================================================================
# WLS state estimation
# =============================================================================
@dataclass
class EstimationResult:
    v: np.ndarray
    converged: bool
    iterations: int
    residual: np.ndarray = field(default_factory=lambda: np.array([]))
    objective: float = float("nan")      # J(x) = sum (r_i / sigma_i)^2
    chi2_threshold: float = float("nan")
    bad_data_detected: bool = False
    largest_normalised_residual: float = float("nan")
    suspect_measurement: int = -1


def estimate(
    ppc: dict[str, Any],
    meas: MeasurementSet,
    tol: float = 1e-6,
    max_iter: int = 30,
    alpha: float = 0.01,
) -> EstimationResult:
    """Weighted least-squares state estimation by Gauss-Newton.

        x_{k+1} = x_k + (H' W H)^-1 H' W (z - h(x_k))

    Then a chi-squared test on J(x): if the weighted residual is larger than the
    distribution predicts for (m - n) degrees of freedom, something is wrong with
    the data -- a failed sensor, a wrong topology, or an attack.
    """
    ref, pv, pq = bus_types(ppc)
    n = ppc["bus"].shape[0]
    non_ref = np.array([i for i in range(n) if i not in set(ref.tolist())])

    v = np.ones(n, dtype=complex)
    v[ref] = ppc["bus"][ref, 7] * np.exp(1j * np.deg2rad(ppc["bus"][ref, 8]))

    w = meas.weights
    converged = False
    it = 0

    for it in range(1, max_iter + 1):
        h_x = measurement_function(ppc, meas, v)
        r = meas.values - h_x
        jac = jacobian(ppc, meas, v)

        gain = jac.T @ (w[:, None] * jac)
        rhs = jac.T @ (w * r)
        try:
            dx = np.linalg.solve(gain + 1e-10 * np.eye(gain.shape[0]), rhs)
        except np.linalg.LinAlgError:
            return EstimationResult(v, False, it)

        va, vm = np.angle(v), np.abs(v)
        va[non_ref] += dx[:non_ref.size]
        vm += dx[non_ref.size:]
        vm = np.clip(vm, 0.3, 2.0)
        v = vm * np.exp(1j * va)

        if np.max(np.abs(dx)) < tol:
            converged = True
            break

    h_x = measurement_function(ppc, meas, v)
    r = meas.values - h_x
    objective = float(np.sum(w * r ** 2))

    dof = max(meas.n - (non_ref.size + n), 1)
    threshold = float(stats.chi2.ppf(1.0 - alpha, dof))

    # normalised residuals identify WHICH measurement is suspect
    jac = jacobian(ppc, meas, v)
    gain = jac.T @ (w[:, None] * jac)
    try:
        omega = np.diag(1.0 / w) - jac @ np.linalg.solve(gain, jac.T)
        denom = np.sqrt(np.clip(np.diag(omega), 1e-12, None))
        r_norm = np.abs(r) / denom
    except np.linalg.LinAlgError:
        r_norm = np.abs(r) / meas.sigma

    return EstimationResult(
        v=v, converged=converged, iterations=it, residual=r,
        objective=objective, chi2_threshold=threshold,
        bad_data_detected=bool(objective > threshold),
        largest_normalised_residual=float(np.max(r_norm)),
        suspect_measurement=int(np.argmax(r_norm)),
    )


# =============================================================================
# false data injection
# =============================================================================
def random_attack(
    meas: MeasurementSet, n_targets: int = 3, magnitude: float = 0.15,
    seed: int = 7,
) -> np.ndarray:
    """Naive attack: perturb a few measurements arbitrarily.

    The baseline a defender should catch. Included so the stealthy attack has
    something to be compared against.
    """
    rng = np.random.default_rng(seed)
    a = np.zeros(meas.n)
    targets = rng.choice(meas.n, size=min(n_targets, meas.n), replace=False)
    a[targets] = magnitude * rng.choice([-1.0, 1.0], size=targets.size)
    return a


def stealthy_attack(
    ppc: dict[str, Any],
    meas: MeasurementSet,
    v: np.ndarray,
    target_bus: int,
    angle_shift_rad: float = 0.05,
    vm_shift: float = 0.0,
    linear: bool = False,
) -> np.ndarray:
    """Construct a residual-preserving false data injection attack.

    The attacker picks a false state change `c` -- here a voltage angle (and
    optionally magnitude) shift at one bus -- and injects the measurement
    perturbation that is exactly consistent with it.

    EXACT (AC) FORM, the default:

        a = h(x + c) - h(x)

    Every measurement is moved to precisely what it would read if the state
    really were x + c. The estimator converges to x + c with an IDENTICAL
    residual, so J(x) is unchanged and the chi-squared test cannot fire. Not
    because the threshold is mistuned -- because there is no residual to test.

    LINEARISED FORM (`linear=True`), a = H c:

        The textbook DC-model construction. It is NOT exactly stealthy against
        an AC estimator: h is nonlinear, so h(x+c) - h(x) - Hc leaves a
        second-order residual that grows with |c|. Measured on IEEE 14-bus with
        a 2.9 deg shift, this raises J from 62 to 110 against a threshold of 82
        -- i.e. the linearised attack IS caught. Kept here because the
        difference between the two is the point: an attacker who models the grid
        properly defeats the defence, and one who does not, does not.

    Detecting the exact form requires information the residual does not contain:
    redundant PMU physics, or cryptographically attested measurements.
    """
    idx = bus_index(ppc)
    i = idx[int(target_bus)]

    ref, pv, pq = bus_types(ppc)
    n = ppc["bus"].shape[0]
    non_ref = np.array([k for k in range(n) if k not in set(ref.tolist())])

    if linear:
        c = np.zeros(non_ref.size + n)
        pos = np.flatnonzero(non_ref == i)
        if pos.size:
            c[pos[0]] = angle_shift_rad
        c[non_ref.size + i] = vm_shift
        return jacobian(ppc, meas, v) @ c

    # exact: move the state, then take the difference of the true measurement
    # function -- no linearisation, so no linearisation residual
    v_attacked = v.copy()
    if i not in set(ref.tolist()):
        v_attacked[i] = (abs(v[i]) + vm_shift) * np.exp(
            1j * (np.angle(v[i]) + angle_shift_rad))
    else:
        v_attacked[i] = (abs(v[i]) + vm_shift) * np.exp(1j * np.angle(v[i]))

    return measurement_function(ppc, meas, v_attacked) - \
        measurement_function(ppc, meas, v)


def apply_attack(meas: MeasurementSet, a: np.ndarray) -> MeasurementSet:
    """Return a copy of the measurement set with the attack vector added."""
    return MeasurementSet(
        types=list(meas.types), indices=list(meas.indices),
        values=meas.values + a, truth=meas.truth.copy(), sigma=meas.sigma.copy(),
    )
