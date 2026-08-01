"""Transient stability: multi-machine swing simulation and critical clearing time.

This closes the one hypothesis the phasor studies could not reach --
`gfm_extends_cct` in config/tolerances.yaml. Everything up to WP6 is a
steady-state or quasi-steady-state question. Whether a grid survives a fault is
a DYNAMIC one, and it is where the difference between grid-following and
grid-forming control shows up most sharply.

Model (classical / "transient" model, the standard first-cut TSA formulation):

  SYNCHRONOUS MACHINE   constant voltage E' behind transient reactance X'd,
                        rotor angle governed by the swing equation

                            d(delta)/dt = omega - omega_s
                        2H/omega_s d(omega)/dt = Pm - Pe - D (omega - omega_s)

  GRID-FORMING (GFM)    also a voltage source behind an impedance, and with
                        virtual-synchronous-machine control it obeys the SAME
                        swing equation with a virtual inertia constant. So a GFM
                        converter participates in the electromechanical dynamics
                        -- that is the entire point of the control.

  GRID-FOLLOWING (GFL)  a current source with NO inertia and no rotor angle. It
                        contributes nothing to the swing dynamics and is folded
                        into the network as a constant-impedance injection.

The network is Kron-reduced to the internal buses of the dynamic sources, which
is what makes the classical model cheap enough to bisect on clearing time.

STATED LIMITATIONS -- read before using a CCT from this module:

  1. Constant E', no exciter, no governor, no saliency, loads as constant
     impedance. A screening model, not an EMT study.

  2. CCT IS NOT COMPARABLE ACROSS SCENARIOS WITH DIFFERENT NUMBERS OF DYNAMIC
     UNITS. This is the important one, and it is easy to get wrong.

     The stability criterion is rotor-angle separation, which is defined only
     over the units that HAVE a rotor angle. Converting a synchronous machine to
     a grid-following converter removes it from the swing model entirely -- so
     the measurement loses the very failure mode it exists to detect. Measured
     on IEEE 39-bus: CCT "improves" from 0.121 s at 0% IBR to 0.977 s at 40.6%
     grid-following penetration. That is an artefact of counting fewer rotors,
     not a physical improvement, and reporting it as one would be wrong.

     A grid-following converter's actual instability mechanism is PLL loss of
     synchronisation -- a converter CONTROL phenomenon on a millisecond
     timescale, which an electromechanical model has no representation of. This
     is a concrete demonstration of why EMT analysis becomes mandatory in
     converter-dominated grids, reached independently of the fault-analysis
     argument in WP3.

     USE `virtual_inertia_sweep` INSTEAD. Holding the unit count and topology
     fixed and varying only the grid-forming virtual inertia is a controlled
     comparison, and it is the one this module can answer honestly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from . import ybus as yb
from .config import scenarios
from .ppc import (
    BR_STATUS, GEN_BUS, GEN_STATUS, PD, PG, QD, bus_index, gen_at_bus,
)


# =============================================================================
# network reduction
# =============================================================================
@dataclass
class DynamicNetwork:
    """Network reduced to the internal buses of the dynamic (voltage-source) units."""

    y_reduced: np.ndarray          # admittance between internal source buses
    e_internal: np.ndarray         # complex internal EMFs (constant magnitude)
    h: np.ndarray                  # inertia constants, seconds, on system base
    p_mech: np.ndarray             # mechanical power, per-unit
    damping: np.ndarray            # damping coefficients
    gen_rows: list[int]
    is_gfm: np.ndarray             # which dynamic units are grid-forming
    delta0: np.ndarray             # initial rotor angles, radians

    @property
    def n(self) -> int:
        return self.y_reduced.shape[0]


def _augmented_ybus(
    ppc: dict[str, Any],
    v: np.ndarray,
    dyn_buses: np.ndarray,
    x_transient: np.ndarray,
    gfl_buses: np.ndarray,
    fault_bus: int | None = None,
    outaged_branch: int | None = None,
) -> np.ndarray:
    """Ybus with loads as constant impedance and dynamic units as internal nodes.

    Loads, and grid-following converters, become fixed shunt admittances derived
    from the pre-fault solution. A GFL unit has no rotor and no internal EMF, so
    representing it as a constant-impedance injection is the consistent choice --
    it can neither absorb nor release rotational energy.
    """
    work = {**ppc, "branch": ppc["branch"].copy()}
    if outaged_branch is not None:
        work["branch"][outaged_branch, BR_STATUS] = 0

    n = ppc["bus"].shape[0]
    y = yb.build(work).astype(complex)
    diag = np.arange(n)

    # loads as constant impedance from the prefault operating point
    s_load = (ppc["bus"][:, PD] + 1j * ppc["bus"][:, QD]) / ppc["baseMVA"]
    vm2 = np.abs(v) ** 2
    vm2[vm2 < 1e-9] = 1.0
    y[diag, diag] += np.conj(s_load) / vm2

    # grid-following injections as constant impedance (negative load)
    gens = gen_at_bus(ppc)
    for i in gfl_buses:
        rows = gens.get(int(i), [])
        if not rows:
            continue
        p = float(ppc["gen"][rows, PG].sum()) / ppc["baseMVA"]
        s_inj = complex(p, 0.0)
        y[i, i] += np.conj(-s_inj) / vm2[i]

    # bolted three-phase fault: very large shunt to ground
    if fault_bus is not None:
        y[fault_bus, fault_bus] += 1e7

    # dynamic units: add an internal node behind X'd
    n_dyn = dyn_buses.size
    total = n + n_dyn
    y_aug = np.zeros((total, total), dtype=complex)
    y_aug[:n, :n] = y

    for k, (i, x) in enumerate(zip(dyn_buses, x_transient)):
        node = n + k
        y_series = 1.0 / complex(0.0, float(x))
        y_aug[node, node] += y_series
        y_aug[i, i] += y_series
        y_aug[node, i] -= y_series
        y_aug[i, node] -= y_series

    return y_aug


def _kron_reduce(y_aug: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """Eliminate every node except `keep` (the internal source nodes)."""
    n = y_aug.shape[0]
    drop = np.array([i for i in range(n) if i not in set(keep.tolist())], dtype=int)
    if drop.size == 0:
        return y_aug[np.ix_(keep, keep)]

    y_kk = y_aug[np.ix_(keep, keep)]
    y_kd = y_aug[np.ix_(keep, drop)]
    y_dk = y_aug[np.ix_(drop, keep)]
    y_dd = y_aug[np.ix_(drop, drop)]
    return y_kk - y_kd @ np.linalg.solve(y_dd, y_dk)


def build_dynamic_network(
    ppc: dict[str, Any],
    v: np.ndarray,
    ibr_gens: np.ndarray | None = None,
    gfm_gens: np.ndarray | None = None,
    h_seconds: float | None = None,
) -> DynamicNetwork:
    """Assemble the classical transient-stability model for a scenario.

    Dynamic (voltage-source) units are the synchronous machines PLUS the
    grid-forming converters. Grid-following converters are excluded: they have
    no rotor angle to swing.
    """
    cfg = scenarios()["dynamics"]
    md = scenarios()["machine_data"]
    h_default = float(cfg.get("h_default_s", 4.0)) if h_seconds is None else h_seconds
    h_virtual = float(cfg["gfm"].get("inertia_const_s", 4.0))
    x_dp = float(md.get("xdp_pu", 0.30))          # transient reactance X'd
    zv = cfg["gfm"]["virtual_impedance_pu"]

    from .ibr import inverter_rating_pu

    # published per-machine data where the case is a standard dynamic benchmark
    dyn_data = scenarios().get("machine_dynamics", {}).get(ppc.get("case", ""), {})
    dyn_data = {int(k): v for k, v in dyn_data.items()}

    idx = bus_index(ppc)
    ibr_set = set(np.asarray(ibr_gens, dtype=int).tolist()) if ibr_gens is not None else set()
    gfm_set = set(np.asarray(gfm_gens, dtype=int).tolist()) if gfm_gens is not None else set()
    gfm_set &= ibr_set
    gfl_set = ibr_set - gfm_set

    live = np.flatnonzero(ppc["gen"][:, GEN_STATUS] > 0)
    dyn_rows = [int(g) for g in live if int(g) not in gfl_set]
    gfl_buses = np.array(sorted({idx[int(ppc["gen"][g, GEN_BUS])] for g in gfl_set}), dtype=int)

    dyn_buses, x_tr, h_list, is_gfm = [], [], [], []
    base = ppc["baseMVA"]
    for g in dyn_rows:
        bus_no = int(ppc["gen"][g, GEN_BUS])
        bus_row = idx[bus_no]
        rating = ppc["gen_mva_base"][g]
        if not np.isfinite(rating) or rating <= 0:
            rating = max(abs(float(ppc["gen"][g, PG])), base)

        published = dyn_data.get(bus_no)

        if g in gfm_set:
            # a grid-forming converter behind its virtual impedance, with the
            # virtual inertia its VSM control synthesises
            s_pu = inverter_rating_pu(ppc, [g])
            x = float(zv["x"]) / s_pu if s_pu > 0 else x_dp
            h_unit = h_virtual * rating / base
            is_gfm.append(True)
        elif published is not None:
            # published benchmark data, already on the system base
            x = float(published["xdp"])
            h_unit = float(published["h"])
            is_gfm.append(False)
        else:
            x = x_dp * base / rating          # machine base -> system base
            h_unit = h_default * rating / base
            is_gfm.append(False)

        dyn_buses.append(bus_row)
        x_tr.append(x)
        h_list.append(h_unit)

    dyn_buses = np.array(dyn_buses, dtype=int)
    x_tr = np.array(x_tr, dtype=float)

    # internal EMFs from the prefault solution: E' = V + jX' I
    ybus = yb.build(ppc)
    s_bus = v * np.conj(ybus @ v)
    e_int = np.zeros(dyn_buses.size, dtype=complex)
    p_mech = np.zeros(dyn_buses.size)
    for k, (g, i) in enumerate(zip(dyn_rows, dyn_buses)):
        p_gen = float(ppc["gen"][g, PG]) / base
        q_gen = float(s_bus[i].imag) + float(ppc["bus"][i, 3]) / base
        i_gen = np.conj(complex(p_gen, q_gen) / v[i]) if abs(v[i]) > 1e-9 else 0j
        e_int[k] = v[i] + 1j * x_tr[k] * i_gen
        p_mech[k] = p_gen

    n = ppc["bus"].shape[0]
    keep = np.arange(n, n + dyn_buses.size)
    y_aug = _augmented_ybus(ppc, v, dyn_buses, x_tr, gfl_buses)
    y_red = _kron_reduce(y_aug, keep)

    return DynamicNetwork(
        y_reduced=y_red, e_internal=e_int, h=np.array(h_list),
        p_mech=p_mech, damping=np.full(dyn_buses.size, 2.0),
        gen_rows=dyn_rows, is_gfm=np.array(is_gfm, dtype=bool),
        delta0=np.angle(e_int),
    )


def reduced_for_condition(
    ppc: dict[str, Any], v: np.ndarray, net: DynamicNetwork,
    ibr_gens: np.ndarray | None, gfm_gens: np.ndarray | None,
    fault_bus: int | None = None, outaged_branch: int | None = None,
) -> np.ndarray:
    """Reduced admittance matrix for the faulted or post-fault topology."""
    idx = bus_index(ppc)
    ibr_set = set(np.asarray(ibr_gens, dtype=int).tolist()) if ibr_gens is not None else set()
    gfm_set = set(np.asarray(gfm_gens, dtype=int).tolist()) if gfm_gens is not None else set()
    gfl_set = ibr_set - (gfm_set & ibr_set)
    gfl_buses = np.array(sorted({idx[int(ppc["gen"][g, GEN_BUS])] for g in gfl_set}), dtype=int)

    dyn_buses = np.array([idx[int(ppc["gen"][g, GEN_BUS])] for g in net.gen_rows], dtype=int)

    cfg = scenarios()["dynamics"]
    md = scenarios()["machine_data"]
    from .ibr import inverter_rating_pu
    zv = cfg["gfm"]["virtual_impedance_pu"]
    x_dp = float(md.get("xdp_pu", 0.30))
    base = ppc["baseMVA"]

    dyn_data = scenarios().get("machine_dynamics", {}).get(ppc.get("case", ""), {})
    dyn_data = {int(k): v for k, v in dyn_data.items()}

    x_tr = []
    for k, g in enumerate(net.gen_rows):
        bus_no = int(ppc["gen"][g, GEN_BUS])
        rating = ppc["gen_mva_base"][g]
        if not np.isfinite(rating) or rating <= 0:
            rating = max(abs(float(ppc["gen"][g, PG])), base)
        published = dyn_data.get(bus_no)
        if net.is_gfm[k]:
            s_pu = inverter_rating_pu(ppc, [g])
            x_tr.append(float(zv["x"]) / s_pu if s_pu > 0 else x_dp)
        elif published is not None:
            x_tr.append(float(published["xdp"]))
        else:
            x_tr.append(x_dp * base / rating)

    n = ppc["bus"].shape[0]
    keep = np.arange(n, n + dyn_buses.size)
    y_aug = _augmented_ybus(ppc, v, dyn_buses, np.array(x_tr), gfl_buses,
                            fault_bus=fault_bus, outaged_branch=outaged_branch)
    return _kron_reduce(y_aug, keep)


# =============================================================================
# swing integration
# =============================================================================
@dataclass
class SwingResult:
    t: np.ndarray
    delta: np.ndarray              # (n_steps, n_units) radians
    omega: np.ndarray              # per-unit deviation
    stable: bool
    max_separation_deg: float
    note: str = ""
    coi_delta: np.ndarray = field(default_factory=lambda: np.array([]))


def _electrical_power(e: np.ndarray, delta: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pe_i = sum_j |Ei||Ej| (G_ij cos(di-dj) + B_ij sin(di-dj))."""
    mag = np.abs(e)
    v = mag * np.exp(1j * delta)
    s = v * np.conj(y @ v)
    return np.real(s)


def simulate(
    net: DynamicNetwork,
    y_fault: np.ndarray,
    y_post: np.ndarray,
    t_clear: float,
    t_end: float = 5.0,
    dt: float = 0.002,
    f0: float = 60.0,
    separation_limit_deg: float = 180.0,
) -> SwingResult:
    """Integrate the swing equations through fault-on and post-fault periods.

    RK4 with a fixed step. Instability is declared when any rotor angle departs
    from the centre-of-inertia by more than `separation_limit_deg` -- the
    standard loss-of-synchronism criterion.
    """
    omega_s = 2.0 * np.pi * f0
    n_steps = int(t_end / dt)
    n = net.n

    delta = net.delta0.copy()
    omega = np.zeros(n)
    mag = np.abs(net.e_internal)

    t_hist = np.zeros(n_steps + 1)
    d_hist = np.zeros((n_steps + 1, n))
    w_hist = np.zeros((n_steps + 1, n))
    d_hist[0] = delta
    w_hist[0] = omega

    h_tot = net.h.sum()
    stable = True
    max_sep = 0.0
    note = ""

    def deriv(d: np.ndarray, w: np.ndarray, y: np.ndarray):
        pe = _electrical_power(mag * np.exp(1j * 0), d, y) if False else \
            _electrical_power(net.e_internal, d, y)
        dd = w * omega_s
        dw = (net.p_mech - pe - net.damping * w) / (2.0 * net.h)
        return dd, dw

    for k in range(n_steps):
        t = k * dt
        y = y_fault if t < t_clear else y_post

        k1d, k1w = deriv(delta, omega, y)
        k2d, k2w = deriv(delta + 0.5 * dt * k1d, omega + 0.5 * dt * k1w, y)
        k3d, k3w = deriv(delta + 0.5 * dt * k2d, omega + 0.5 * dt * k2w, y)
        k4d, k4w = deriv(delta + dt * k3d, omega + dt * k3w, y)

        delta = delta + dt / 6.0 * (k1d + 2 * k2d + 2 * k3d + k4d)
        omega = omega + dt / 6.0 * (k1w + 2 * k2w + 2 * k3w + k4w)

        t_hist[k + 1] = t + dt
        d_hist[k + 1] = delta
        w_hist[k + 1] = omega

        coi = float(np.sum(net.h * delta) / h_tot) if h_tot > 0 else 0.0
        sep = float(np.max(np.abs(np.degrees(delta - coi))))
        max_sep = max(max_sep, sep)

        if not np.isfinite(delta).all():
            stable, note = False, "numerical divergence"
            break
        if sep > separation_limit_deg:
            stable, note = False, f"loss of synchronism at t={t + dt:.3f}s"
            break

    used = k + 2 if not stable else n_steps + 1
    coi_series = (d_hist[:used] @ net.h) / h_tot if h_tot > 0 else np.zeros(used)
    return SwingResult(
        t=t_hist[:used], delta=d_hist[:used], omega=w_hist[:used],
        stable=stable, max_separation_deg=max_sep, note=note,
        coi_delta=coi_series,
    )


def critical_clearing_time(
    ppc: dict[str, Any],
    v: np.ndarray,
    net: DynamicNetwork,
    fault_bus: int,
    ibr_gens: np.ndarray | None = None,
    gfm_gens: np.ndarray | None = None,
    outaged_branch: int | None = None,
    t_min: float = 0.0,
    t_max: float = 0.6,
    tol: float = 0.005,
    t_end: float = 5.0,
) -> dict[str, Any]:
    """Bisection search for the longest clearing time the system survives.

    CCT is the standard measure of transient stability margin, and the quantity
    the GFM-vs-GFL comparison turns on: more effective inertia and a stiffer
    voltage source behind the converter should let the system tolerate a longer
    fault before losing synchronism.
    """
    idx = bus_index(ppc)
    fb = idx[int(fault_bus)]

    y_fault = reduced_for_condition(ppc, v, net, ibr_gens, gfm_gens, fault_bus=fb)
    y_post = reduced_for_condition(ppc, v, net, ibr_gens, gfm_gens,
                                   outaged_branch=outaged_branch)

    # a system unstable even with instantaneous clearing has no CCT
    probe = simulate(net, y_fault, y_post, t_clear=t_min, t_end=t_end)
    if not probe.stable:
        return {"cct_s": 0.0, "bracketed": False,
                "note": "unstable even at zero clearing time"}

    probe = simulate(net, y_fault, y_post, t_clear=t_max, t_end=t_end)
    if probe.stable:
        return {"cct_s": t_max, "bracketed": False,
                "note": f"still stable at t_max={t_max}s"}

    lo, hi = t_min, t_max
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if simulate(net, y_fault, y_post, t_clear=mid, t_end=t_end).stable:
            lo = mid
        else:
            hi = mid

    return {"cct_s": lo, "bracketed": True, "note": "",
            "n_dynamic_units": net.n, "h_total_s": float(net.h.sum())}


def virtual_inertia_sweep(
    ppc: dict[str, Any],
    v: np.ndarray,
    fault_bus: int,
    ibr_gens: np.ndarray,
    gfm_gens: np.ndarray,
    h_virtual_values: list[float],
    outaged_branch: int | None = None,
    t_max: float = 1.0,
) -> list[dict[str, Any]]:
    """CCT vs grid-forming virtual inertia, at FIXED topology and unit count.

    This is the controlled version of the GFM question. The set of dynamic units,
    the network, the fault and the dispatch are all held constant; the ONLY thing
    that changes is how much virtual inertia the grid-forming converters
    synthesise. Any change in CCT is therefore attributable to synthetic inertia
    and nothing else -- unlike a GFM-vs-GFL comparison, which silently changes
    how many rotor angles exist and so changes what the metric even measures.
    """
    rows = []
    base_net = build_dynamic_network(ppc, v, ibr_gens=ibr_gens, gfm_gens=gfm_gens)
    gfm_mask = base_net.is_gfm

    for h_v in h_virtual_values:
        net = build_dynamic_network(ppc, v, ibr_gens=ibr_gens, gfm_gens=gfm_gens)
        net.damping[:] = 0.0
        # scale only the grid-forming units' inertia
        if gfm_mask.any():
            base_h = build_dynamic_network(
                ppc, v, ibr_gens=ibr_gens, gfm_gens=gfm_gens).h[gfm_mask]
            h_nominal = float(scenarios()["dynamics"]["gfm"]["inertia_const_s"])
            net.h[gfm_mask] = base_h * (h_v / h_nominal) if h_nominal > 0 else base_h

        res = critical_clearing_time(
            ppc, v, net, fault_bus, ibr_gens=ibr_gens, gfm_gens=gfm_gens,
            outaged_branch=outaged_branch, t_max=t_max,
        )
        rows.append({
            "h_virtual_s": h_v,
            "cct_s": res["cct_s"],
            "bracketed": res["bracketed"],
            "n_dynamic_units": net.n,
            "h_total_s": float(net.h.sum()),
            "h_gfm_total_s": float(net.h[gfm_mask].sum()) if gfm_mask.any() else 0.0,
            "note": res.get("note", ""),
        })
    return rows
