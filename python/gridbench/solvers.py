"""Power flow solvers written from first principles: Newton-Raphson, Gauss-Seidel,
and Fast-Decoupled (XB / BX).

These exist to be *compared*, not just to converge. Every solver returns the same
`PFResult` carrying iteration count, the full mismatch history, wall-clock time
and (for NR) the Jacobian condition number -- which is the quantity that actually
explains why NR degrades on ill-conditioned, low-inertia, high-IBR cases.

Reference implementations (MATPOWER `runpf`, pandapower `runpp`) are used only as
an oracle in the test suite; nothing here calls them.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from . import ybus as yb
from .ppc import (
    BR_B, BR_R, BS, GEN_BUS, GEN_STATUS, PQ, PV, QG, QMAX, QMIN, REF, SHIFT,
    TAP, VG, VA, VM, bus_index, bus_injections, bus_types, gen_at_bus,
)

Algorithm = Literal["nr", "gs", "fdxb", "fdbx"]


@dataclass
class PFResult:
    """Outcome of one power flow solve."""

    algorithm: str
    converged: bool
    iterations: int
    v: np.ndarray                       # complex bus voltages, per-unit
    mismatch_history: list[float] = field(default_factory=list)
    elapsed_s: float = 0.0
    jacobian_cond: float = float("nan")  # NR only; condition number at solution
    q_limits_hit: list[int] = field(default_factory=list)  # buses switched PV->PQ
    note: str = ""

    @property
    def vm(self) -> np.ndarray:
        return np.abs(self.v)

    @property
    def va_deg(self) -> np.ndarray:
        return np.degrees(np.angle(self.v))

    @property
    def final_mismatch(self) -> float:
        return self.mismatch_history[-1] if self.mismatch_history else float("nan")

    def summary(self) -> str:
        state = "converged" if self.converged else "DIVERGED"
        return (
            f"{self.algorithm:<6} {state:>9}  it={self.iterations:>4}  "
            f"mismatch={self.final_mismatch:.3e}  t={self.elapsed_s * 1e3:7.2f} ms"
        )


# -----------------------------------------------------------------------------
# shared helpers
# -----------------------------------------------------------------------------
def _flat_start(ppc: dict[str, Any]) -> np.ndarray:
    """Flat start with generator setpoints applied at PV/slack buses."""
    bus = ppc["bus"]
    v = bus[:, VM] * np.exp(1j * np.deg2rad(bus[:, VA]))
    idx = bus_index(ppc)
    for row in ppc["gen"]:
        if row[GEN_STATUS] > 0:
            i = idx[int(row[GEN_BUS])]
            v[i] = row[VG] * np.exp(1j * np.angle(v[i]))
    return v


def _power_mismatch(
    v: np.ndarray, ybus: np.ndarray, sbus: np.ndarray,
    pvpq: np.ndarray, pq: np.ndarray,
) -> np.ndarray:
    """Real mismatch vector [dP(pv,pq); dQ(pq)]."""
    mis = v * np.conj(ybus @ v) - sbus
    return np.concatenate([mis[pvpq].real, mis[pq].imag])


def _dsbus_dv(ybus: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Partial derivatives of bus injections w.r.t. voltage magnitude and angle.

    Standard polar-coordinate result:
        dS/d|V| = diag(V/|V|) conj(diag(Ibus)) + diag(V) conj(Ybus diag(V/|V|))
        dS/dtheta = j diag(V) conj(diag(Ibus) - Ybus diag(V))
    """
    ibus = ybus @ v
    vnorm = v / np.abs(v)
    ds_dvm = np.diag(vnorm) @ np.conj(np.diag(ibus)) + np.diag(v) @ np.conj(ybus @ np.diag(vnorm))
    ds_dva = 1j * np.diag(v) @ np.conj(np.diag(ibus) - ybus @ np.diag(v))
    return ds_dvm, ds_dva


def _build_jacobian(
    ybus: np.ndarray, v: np.ndarray, pvpq: np.ndarray, pq: np.ndarray,
) -> np.ndarray:
    """Full polar Newton-Raphson Jacobian [[dP/dth, dP/dVm], [dQ/dth, dQ/dVm]]."""
    ds_dvm, ds_dva = _dsbus_dv(ybus, v)
    j11 = ds_dva[np.ix_(pvpq, pvpq)].real
    j12 = ds_dvm[np.ix_(pvpq, pq)].real
    j21 = ds_dva[np.ix_(pq, pvpq)].imag
    j22 = ds_dvm[np.ix_(pq, pq)].imag
    return np.block([[j11, j12], [j21, j22]])


def _gen_q_injection(
    ppc: dict[str, Any], v: np.ndarray, ybus: np.ndarray, i: int,
) -> float:
    """Reactive power (per-unit) the generators at internal bus `i` must supply."""
    q_inj = float((v * np.conj(ybus @ v))[i].imag)      # net injection at the bus
    q_load = float(ppc["bus"][i, 3] / ppc["baseMVA"])   # QD column
    return q_inj + q_load


# -----------------------------------------------------------------------------
# Newton-Raphson
# -----------------------------------------------------------------------------
def newton_raphson(
    ppc: dict[str, Any],
    tol: float = 1e-10,
    max_iter: int = 50,
    enforce_q_limits: bool = False,
    ybus: np.ndarray | None = None,
    v0: np.ndarray | None = None,
) -> PFResult:
    """Full Newton-Raphson power flow in polar coordinates.

    Quadratic convergence: on a well-conditioned IEEE case this reaches 1e-10
    in 3-5 iterations regardless of system size. With `enforce_q_limits`, PV
    buses whose reactive output exceeds Qmax/Qmin are converted to PQ and the
    solve restarts -- this matters for WP2, because the whole point of an IBR is
    that its reactive capability is bounded by a converter current limit.
    """
    t0 = time.perf_counter()
    ybus = yb.build(ppc) if ybus is None else ybus
    sbus = bus_injections(ppc)
    ref, pv, pq = bus_types(ppc)
    # A warm start matters at high IBR penetration: with much of the voltage
    # support removed, the solution can sit far from flat, and NR's basin of
    # attraction shrinks. Starting from the previous scenario's answer separates
    # genuine non-existence of a solution from mere initialisation failure.
    v = _flat_start(ppc) if v0 is None else np.asarray(v0, dtype=complex).copy()

    switched: list[int] = []
    history: list[float] = []
    total_iter = 0
    converged = False

    for _outer in range(10 if enforce_q_limits else 1):
        pvpq = np.sort(np.concatenate([pv, pq]))
        f = _power_mismatch(v, ybus, sbus, pvpq, pq)
        history.append(float(np.max(np.abs(f))) if f.size else 0.0)
        converged = bool(f.size == 0 or np.max(np.abs(f)) < tol)

        it = 0
        while not converged and it < max_iter:
            it += 1
            total_iter += 1
            jac = _build_jacobian(ybus, v, pvpq, pq)
            try:
                dx = np.linalg.solve(jac, -f)
            except np.linalg.LinAlgError:
                return PFResult("nr", False, total_iter, v, history,
                                time.perf_counter() - t0,
                                note="singular Jacobian")
            va, vm = np.angle(v), np.abs(v)
            va[pvpq] += dx[: pvpq.size]
            vm[pq] += dx[pvpq.size:]
            v = vm * np.exp(1j * va)

            f = _power_mismatch(v, ybus, sbus, pvpq, pq)
            history.append(float(np.max(np.abs(f))))
            converged = bool(np.max(np.abs(f)) < tol)

        if not enforce_q_limits or not converged:
            break

        # --- Q-limit enforcement: convert violating PV buses to PQ -----------
        violated = False
        for i, gens in gen_at_bus(ppc).items():
            if i not in pv:
                continue
            q = _gen_q_injection(ppc, v, ybus, i) * ppc["baseMVA"]
            qmax = float(ppc["gen"][gens, QMAX].sum())
            qmin = float(ppc["gen"][gens, QMIN].sum())
            if q > qmax + 1e-6 or q < qmin - 1e-6:
                q_fix = min(max(q, qmin), qmax)
                sbus[i] += 1j * (q_fix - q) / ppc["baseMVA"]
                pv = pv[pv != i]
                pq = np.sort(np.append(pq, i))
                switched.append(int(ppc["bus"][i, 0]))
                violated = True
        if not violated:
            break

    cond = float("nan")
    if converged:
        pvpq = np.sort(np.concatenate([pv, pq]))
        try:
            cond = float(np.linalg.cond(_build_jacobian(ybus, v, pvpq, pq)))
        except np.linalg.LinAlgError:
            pass

    return PFResult("nr", converged, total_iter, v, history,
                    time.perf_counter() - t0, cond, switched)


# -----------------------------------------------------------------------------
# Gauss-Seidel
# -----------------------------------------------------------------------------
def gauss_seidel(
    ppc: dict[str, Any],
    tol: float = 1e-8,
    max_iter: int = 5000,
    accel: float = 1.6,
    ybus: np.ndarray | None = None,
) -> PFResult:
    """Gauss-Seidel power flow with an acceleration factor.

    Linear convergence, and the iteration count scales badly with system size --
    that degradation IS the WP1 result, so `max_iter` is set high enough to let
    it finish rather than being cut off early.
    """
    t0 = time.perf_counter()
    ybus = yb.build(ppc) if ybus is None else ybus
    sbus = bus_injections(ppc)
    ref, pv, pq = bus_types(ppc)
    v = _flat_start(ppc)

    n = v.size
    vg = {i: abs(v[i]) for i in pv}
    qmax_pu, qmin_pu = {}, {}
    for i, gens in gen_at_bus(ppc).items():
        qmax_pu[i] = float(ppc["gen"][gens, QMAX].sum()) / ppc["baseMVA"]
        qmin_pu[i] = float(ppc["gen"][gens, QMIN].sum()) / ppc["baseMVA"]

    ref_set = set(ref.tolist())
    history: list[float] = []
    converged = False
    it = 0

    while it < max_iter:
        it += 1
        v_prev = v.copy()
        for i in range(n):
            if i in ref_set:
                continue
            y_row_sum = ybus[i, :] @ v - ybus[i, i] * v[i]

            if i in pv:
                # PV: recompute Q from the present voltage, then respect limits
                q_calc = -np.imag(np.conj(v[i]) * (y_row_sum + ybus[i, i] * v[i]))
                q_load = ppc["bus"][i, 3] / ppc["baseMVA"]
                q_gen = q_calc + q_load
                q_hi, q_lo = qmax_pu.get(i, np.inf), qmin_pu.get(i, -np.inf)
                if q_gen > q_hi:
                    s_i = complex(sbus[i].real, q_hi - q_load)
                elif q_gen < q_lo:
                    s_i = complex(sbus[i].real, q_lo - q_load)
                else:
                    s_i = complex(sbus[i].real, q_calc)
                    v_new = (np.conj(s_i) / np.conj(v[i]) - y_row_sum) / ybus[i, i]
                    # hold the scheduled magnitude, take only the new angle
                    v[i] = vg[i] * np.exp(1j * np.angle(v_new))
                    continue
            else:
                s_i = sbus[i]

            v_new = (np.conj(s_i) / np.conj(v[i]) - y_row_sum) / ybus[i, i]
            v[i] = v[i] + accel * (v_new - v[i])

        delta = float(np.max(np.abs(v - v_prev)))
        history.append(delta)
        if delta < tol:
            converged = True
            break

    return PFResult("gs", converged, it, v, history, time.perf_counter() - t0,
                    note=f"accel={accel}")


# -----------------------------------------------------------------------------
# Fast-decoupled (XB / BX)
# -----------------------------------------------------------------------------
def _make_b_matrices(ppc: dict[str, Any], variant: str) -> tuple[np.ndarray, np.ndarray]:
    """B' and B'' for the fast-decoupled solvers (MATPOWER `makeB` convention)."""
    # --- B': drop line charging, taps, and (XB only) resistance -------------
    p1 = {**ppc, "bus": ppc["bus"].copy(), "branch": ppc["branch"].copy()}
    p1["branch"][:, BR_B] = 0.0
    p1["branch"][:, TAP] = 1.0
    p1["bus"][:, BS] = 0.0
    if variant == "fdxb":
        p1["branch"][:, BR_R] = 0.0
    bp = -np.imag(yb.build(p1))

    # --- B'': drop phase shifters and (BX only) resistance ------------------
    p2 = {**ppc, "bus": ppc["bus"].copy(), "branch": ppc["branch"].copy()}
    p2["branch"][:, SHIFT] = 0.0
    if variant == "fdbx":
        p2["branch"][:, BR_R] = 0.0
    bpp = -np.imag(yb.build(p2))
    return bp, bpp


def fast_decoupled(
    ppc: dict[str, Any],
    variant: Literal["fdxb", "fdbx"] = "fdxb",
    tol: float = 1e-10,
    max_iter: int = 100,
    ybus: np.ndarray | None = None,
) -> PFResult:
    """Fast-decoupled load flow, XB or BX formulation.

    Exploits the weak P-V / Q-theta coupling of transmission networks (high X/R).
    B' and B'' are constant, so they are factorised once -- each iteration is far
    cheaper than NR, at the cost of linear rather than quadratic convergence.

    Worth noting for the report: the high-X/R assumption is what makes this work,
    and it is exactly the assumption that fails on distribution feeders and on
    converter-dominated networks with virtual impedance.
    """
    t0 = time.perf_counter()
    ybus = yb.build(ppc) if ybus is None else ybus
    sbus = bus_injections(ppc)
    ref, pv, pq = bus_types(ppc)
    pvpq = np.sort(np.concatenate([pv, pq]))
    v = _flat_start(ppc)

    bp, bpp = _make_b_matrices(ppc, variant)
    bp_r = bp[np.ix_(pvpq, pvpq)]
    bpp_r = bpp[np.ix_(pq, pq)]

    history: list[float] = []
    converged = False
    it = 0

    while it < max_iter:
        it += 1
        # --- P-theta half ---------------------------------------------------
        mis = (v * np.conj(ybus @ v) - sbus) / np.abs(v)
        dp = mis[pvpq].real
        if dp.size:
            v = np.abs(v) * np.exp(1j * (np.angle(v) + _scatter(-np.linalg.solve(bp_r, dp), pvpq, v.size)))

        # --- Q-V half -------------------------------------------------------
        mis = (v * np.conj(ybus @ v) - sbus) / np.abs(v)
        dq = mis[pq].imag
        if dq.size:
            v = (np.abs(v) + _scatter(-np.linalg.solve(bpp_r, dq), pq, v.size)) * np.exp(1j * np.angle(v))

        f = _power_mismatch(v, ybus, sbus, pvpq, pq)
        history.append(float(np.max(np.abs(f))) if f.size else 0.0)
        if not f.size or np.max(np.abs(f)) < tol:
            converged = True
            break

    return PFResult(variant, converged, it, v, history, time.perf_counter() - t0)


def _scatter(values: np.ndarray, idx: np.ndarray, n: int) -> np.ndarray:
    out = np.zeros(n)
    out[idx] = values
    return out


# -----------------------------------------------------------------------------
# dispatcher
# -----------------------------------------------------------------------------
SOLVERS = {
    "nr": newton_raphson,
    "gs": gauss_seidel,
    "fdxb": lambda ppc, **kw: fast_decoupled(ppc, "fdxb", **kw),
    "fdbx": lambda ppc, **kw: fast_decoupled(ppc, "fdbx", **kw),
}


def solve(ppc: dict[str, Any], algorithm: Algorithm = "nr", **kwargs: Any) -> PFResult:
    """Run one of the hand-written solvers by name."""
    if algorithm not in SOLVERS:
        raise KeyError(f"unknown algorithm {algorithm!r}; known: {sorted(SOLVERS)}")
    return SOLVERS[algorithm](ppc, **kwargs)
