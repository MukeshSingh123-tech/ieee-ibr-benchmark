"""Inverter-based resource modelling: penetration scenarios and IBR-aware power flow.

The classical load flow treats a generator as a PV bus: it holds voltage using
whatever reactive power it likes, bounded only by a fixed Qmax. An inverter does
not behave that way. Its reactive capability is bounded by a CONVERTER CURRENT
LIMIT, which means Qmax is not a constant -- it shrinks as terminal voltage falls
and as active power rises:

        |S| = |V| * |I| <= |V| * Ilim
        Qmax(V, P) = sqrt( (|V| * Ilim)^2 - P^2 )

That single coupling is what makes an IBR power flow nonlinear in a way the
classical one is not, and it is the mechanism behind the WP2 headline result:
the classical model silently permits reactive support the hardware cannot
deliver, and the error grows as the grid gets weaker.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from . import solvers, ybus as yb
from .config import scenarios
from .ppc import (
    GEN_BUS, GEN_STATUS, PG, PQ, PV, QG, QMAX, QMIN, REF, VG, bus_index,
)


# =============================================================================
# penetration scenarios
# =============================================================================
Basis = Literal["capacity_mva", "dispatched_mw"]


def _gen_weight(ppc: dict[str, Any], basis: Basis) -> np.ndarray:
    """Per-generator quantity that penetration is measured in.

    `capacity_mva` (default) matches how interconnection studies state IBR
    penetration -- installed nameplate. It also handles synchronous condensers
    correctly: a machine dispatching 0 MW still contributes fault current and
    voltage support, so replacing it with an inverter is a real change that a
    dispatched-MW basis would score as zero.
    """
    if basis == "dispatched_mw":
        return np.abs(ppc["gen"][:, PG].astype(float))
    w = ppc["gen_mva_base"].astype(float).copy()
    slack = ppc["is_slack_gen"]
    # slack has no nameplate; use its dispatch as a stand-in so it can be ranked
    w[slack] = np.abs(ppc["gen"][slack, PG])
    return np.where(np.isfinite(w) & (w > 0), w, 0.0)


def select_ibr_gens(
    ppc: dict[str, Any],
    penetration_pct: float,
    order: Literal["largest_first", "smallest_first", "random"] = "largest_first",
    preserve_slack: bool = True,
    basis: Basis = "capacity_mva",
    seed: int = 0,
) -> np.ndarray:
    """Choose which generator ROWS are displaced by IBRs at a given penetration.

    Returns the smallest set, in the configured order, whose combined weight
    reaches the target. Machines with zero weight are never selected -- they
    would not advance the target, and greedily absorbing them was a real bug:
    it silently converted every machine in the system at 20% penetration.

    The slack is preserved below 100% because a phasor load flow needs an angle
    reference. At 100% it is converted too but REMAINS the reference bus, which
    is the correct representation of a grid-forming inverter: it is an IBR, and
    it sets the angle. That is the one case where "100% IBR" is representable in
    a phasor tool at all.
    """
    gen = ppc["gen"]
    live = np.flatnonzero(gen[:, GEN_STATUS] > 0)
    slack = ppc["is_slack_gen"]
    weight = _gen_weight(ppc, basis)

    total = float(weight[live].sum())
    if total <= 0 or penetration_pct <= 0:
        return np.array([], dtype=int)

    keep_slack = preserve_slack and penetration_pct < 100.0
    pool = np.array(
        [g for g in live if weight[g] > 0 and not (keep_slack and slack[g])],
        dtype=int,
    )
    if pool.size == 0:
        return np.array([], dtype=int)

    if order == "largest_first":
        pool = pool[np.argsort(-weight[pool])]
    elif order == "smallest_first":
        pool = pool[np.argsort(weight[pool])]
    else:
        pool = np.random.default_rng(seed).permutation(pool)

    target = penetration_pct / 100.0 * total
    chosen, accumulated = [], 0.0
    for g in pool:
        if accumulated >= target - 1e-9:
            break
        chosen.append(int(g))
        accumulated += float(weight[g])
    return np.array(sorted(chosen), dtype=int)


def actual_penetration_pct(
    ppc: dict[str, Any], ibr_gens: np.ndarray, basis: Basis = "capacity_mva",
) -> float:
    """Realised penetration, which is granular because machines are discrete."""
    weight = _gen_weight(ppc, basis)
    live = ppc["gen"][:, GEN_STATUS] > 0
    total = float(weight[live].sum())
    if total <= 0:
        return 0.0
    return float(weight[np.asarray(ibr_gens, dtype=int)].sum()) / total * 100.0


def max_penetration_pct(
    ppc: dict[str, Any], preserve_slack: bool = True, basis: Basis = "capacity_mva",
) -> float:
    """Highest penetration reachable on this case without displacing the slack.

    Worth reporting alongside any sweep: on IEEE 14-bus with a dispatched-MW
    basis this is only 14.7%, because the slack carries 232 of the 272 MW. A
    sweep that silently clips at its ceiling would otherwise look like a result.
    """
    weight = _gen_weight(ppc, basis)
    live = ppc["gen"][:, GEN_STATUS] > 0
    total = float(weight[live].sum())
    if total <= 0:
        return 0.0
    slack = ppc["is_slack_gen"]
    avail = float(weight[live & ~slack].sum()) if preserve_slack else total
    return avail / total * 100.0


def split_gfl_gfm(
    ppc: dict[str, Any],
    ibr_gens: np.ndarray,
    gfm_share_pct: float,
    order: Literal["largest_first", "smallest_first"] = "largest_first",
) -> tuple[np.ndarray, np.ndarray]:
    """Split an IBR fleet into grid-following and grid-forming units.

    Returns (gfl_gens, gfm_gens). Share is measured in MVA capacity, matching
    how a system operator would specify a grid-forming requirement.

    Largest-first by default: if you are going to mandate grid-forming
    capability on part of a fleet, the large plants are where it buys the most
    system strength per unit of cost.
    """
    ibr_gens = np.asarray(ibr_gens, dtype=int)
    if ibr_gens.size == 0 or gfm_share_pct <= 0:
        return ibr_gens, np.array([], dtype=int)

    weight = _gen_weight(ppc, "capacity_mva")
    total = float(weight[ibr_gens].sum())
    if total <= 0:
        return ibr_gens, np.array([], dtype=int)
    if gfm_share_pct >= 100:
        return np.array([], dtype=int), ibr_gens

    pool = ibr_gens[np.argsort(-weight[ibr_gens])] if order == "largest_first" \
        else ibr_gens[np.argsort(weight[ibr_gens])]

    target = gfm_share_pct / 100.0 * total
    gfm, acc = [], 0.0
    for g in pool:
        if acc >= target - 1e-9:
            break
        gfm.append(int(g))
        acc += float(weight[g])

    gfm_arr = np.array(sorted(gfm), dtype=int)
    gfl_arr = np.array(sorted(set(ibr_gens.tolist()) - set(gfm)), dtype=int)
    return gfl_arr, gfm_arr


def add_synchronous_condensers(
    ppc: dict[str, Any],
    buses: list[int],
    mva_each: float,
) -> dict[str, Any]:
    """Install synchronous condensers at the given buses.

    A synchronous condenser is a machine with no prime mover: P = 0, but it
    still contributes inertia, short-circuit current and a voltage reference.
    This is the mitigation utilities actually buy when system strength runs out
    (National Grid, EirGrid and AEMO have all procured them), and it is one of
    the two levers this project tests -- the other being grid-forming control.

    Returns a NEW ppc; the input is not modified.
    """
    out = {**ppc, "gen": ppc["gen"].copy(), "bus": ppc["bus"].copy()}
    idx = bus_index(out)

    new_rows = []
    for b in buses:
        if int(b) not in idx:
            raise KeyError(f"bus {b} not in case")
        row = np.zeros(out["gen"].shape[1])
        row[GEN_BUS] = int(b)
        row[PG] = 0.0                       # no prime mover
        row[QG] = 0.0
        row[QMAX] = mva_each
        row[QMIN] = -mva_each
        row[VG] = float(out["bus"][idx[int(b)], 7]) or 1.0
        row[6] = mva_each                   # MBASE
        row[GEN_STATUS] = 1
        row[8] = 0.0                        # PMAX = 0
        row[9] = 0.0                        # PMIN
        new_rows.append(row)

    if new_rows:
        out["gen"] = np.vstack([out["gen"], np.array(new_rows)])
        out["is_slack_gen"] = np.concatenate(
            [ppc["is_slack_gen"], np.zeros(len(new_rows), dtype=bool)])
        out["gen_mva_base"] = np.concatenate(
            [ppc["gen_mva_base"], np.full(len(new_rows), float(mva_each))])
        # condenser buses must be able to regulate voltage
        for b in buses:
            i = idx[int(b)]
            if out["bus"][i, 1] == PQ:
                out["bus"][i, 1] = PV
    return out


def ibr_buses_from_gens(ppc: dict[str, Any], ibr_gens: np.ndarray) -> np.ndarray:
    """Internal bus indices hosting at least one IBR."""
    idx = bus_index(ppc)
    return np.array(
        sorted({idx[int(b)] for b in ppc["gen"][np.asarray(ibr_gens, dtype=int), GEN_BUS]}),
        dtype=int,
    )


# =============================================================================
# converter capability
# =============================================================================
def volt_var_q_pu(vm: float, s_rating_pu: float) -> float:
    """IEEE 1547-2018 Category B volt-var droop, returning Q in per-unit.

    Positive = injecting reactive power (supporting a low voltage).
    """
    vv = scenarios()["ibr"]["volt_var"]
    if not vv.get("enabled", True):
        return 0.0
    v1, v2, v3, v4 = vv["v1"], vv["v2"], vv["v3"], vv["v4"]
    q1, q4 = vv["q1_pu"], vv["q4_pu"]

    if vm <= v1:
        q = q1
    elif vm < v2:
        q = q1 * (v2 - vm) / (v2 - v1)
    elif vm <= v3:
        q = 0.0
    elif vm < v4:
        q = q4 * (vm - v3) / (v4 - v3)
    else:
        q = q4
    return q * s_rating_pu


def droop_q_pu(vm: float, v_sched: float, s_rating_pu: float) -> float:
    """Voltage-droop reactive command referenced to a SCHEDULED voltage (IEEE 2800).

        Q = -(V - V_sched) / droop * Q_rating

    This is what a transmission-connected plant controller actually does, and it
    has the property the IEEE 1547 curve lacks here: at the scheduled operating
    point the command is zero, so replacing a synchronous machine with an IBR is
    a smooth change rather than an instantaneous demand for hundreds of MVAr.
    """
    vc = scenarios()["ibr"]["volt_control"]
    err = vm - v_sched
    deadband = float(vc.get("deadband_pu", 0.0))
    if abs(err) <= deadband:
        return 0.0
    err -= np.sign(err) * deadband

    droop = float(vc.get("droop_pct", 4.0)) / 100.0
    q_rating = s_rating_pu * float(vc.get("q_rating_fraction", 0.33))
    return float(-err / droop * q_rating)


def reactive_command_pu(vm: float, v_sched: float, s_rating_pu: float) -> float:
    """Steady-state reactive command, dispatched by the configured control mode."""
    mode = scenarios()["ibr"]["volt_control"].get("mode", "droop_to_setpoint")
    if mode == "ieee1547_voltvar":
        return volt_var_q_pu(vm, s_rating_pu)
    return droop_q_pu(vm, v_sched, s_rating_pu)


def q_capability_pu(
    vm: float, p_pu: float, i_limit: float, s_rating_pu: float,
    apply_plant_cap: bool = True,
) -> float:
    """Reactive capability actually available at an IBR terminal, per-unit.

    Two limits apply and the binding one wins:

      converter current   Qmax = sqrt( (V * Ilim * S)^2 - P^2 )
      plant capability    Qmax = q_capability_fraction * S      (~0.95 pf)

    The first is the physics the classical PV bus lacks; the second is the
    interconnection agreement. In steady state the plant cap usually binds
    first, which is why ignoring it overstates available reactive support. As
    voltage collapses during a fault the current-limit term shrinks toward zero
    and takes over -- which is exactly the WP3 mechanism.
    """
    s_max = vm * i_limit * s_rating_pu
    q_current_limit = float(np.sqrt(max(0.0, s_max ** 2 - p_pu ** 2)))

    cfg = scenarios()["ibr"]
    if cfg.get("sizing", "match_machine_capability") == "match_machine_capability":
        # rating already encodes the machine's capability; a second cap would
        # shrink it below the machine it replaces and break the control
        return q_current_limit
    if not apply_plant_cap:
        return q_current_limit
    frac = float(cfg.get("q_capability_fraction", 0.33))
    return min(q_current_limit, frac * s_rating_pu)


def inverter_rating_pu(ppc: dict[str, Any], gen_rows: list[int]) -> float:
    """MVA rating of the inverter replacing the machines in `gen_rows`, per-unit.

    Under `match_machine_capability`, S = sqrt(PMAX^2 + QMAX^2). Substituting
    that into Qmax(V, P) at V = 1.0, Ilim = 1.0, P = PMAX returns exactly QMAX --
    so the inverter starts life with precisely the reactive capability of the
    machine it displaces, and every subsequent difference is attributable to the
    voltage and power dependence of the converter limit.
    """
    from .ppc import PMAX, QMAX as QMAX_COL

    base = ppc["baseMVA"]
    cfg = scenarios()["ibr"]
    gen = ppc["gen"]

    if cfg.get("sizing", "match_machine_capability") == "match_machine_capability":
        pmax = gen[gen_rows, PMAX].astype(float)
        qmax = gen[gen_rows, QMAX_COL].astype(float)
        # the slack's PMAX is a sentinel (1e9); fall back to its dispatch
        pmax = np.where(np.isfinite(pmax) & (pmax < 1e6), pmax,
                        np.abs(gen[gen_rows, PG]))
        qmax = np.where(np.isfinite(qmax) & (np.abs(qmax) < 1e6), np.abs(qmax), 0.0)
        rating = float(np.sqrt((pmax ** 2 + qmax ** 2).sum()))
    else:
        rating = float(np.nansum(ppc["gen_mva_base"][gen_rows]))

    if not np.isfinite(rating) or rating <= 0:
        rating = max(float(np.abs(gen[gen_rows, PG]).sum()), 1.0)
    return rating / base


def scale_load(ppc: dict[str, Any], factor: float) -> dict[str, Any]:
    """Return a copy of the case with all bus load scaled by `factor`.

    Generation is NOT rescaled: the slack absorbs the imbalance, which is the
    standard way of stressing a case toward its voltage-stability limit.
    """
    out = {**ppc, "bus": ppc["bus"].copy()}
    out["bus"][:, 2] *= factor      # PD
    out["bus"][:, 3] *= factor      # QD
    return out


# =============================================================================
# IBR-aware power flow
# =============================================================================
@dataclass
class IBRPowerFlowResult:
    converged: bool
    outer_iterations: int
    inner_iterations: int
    v: np.ndarray
    q_ibr_mvar: dict[int, float] = field(default_factory=dict)
    q_classical_mvar: dict[int, float] = field(default_factory=dict)
    limited_buses: list[int] = field(default_factory=list)
    penetration_pct: float = 0.0
    note: str = ""

    @property
    def vm(self) -> np.ndarray:
        return np.abs(self.v)

    @property
    def va_deg(self) -> np.ndarray:
        return np.degrees(np.angle(self.v))


def ibr_powerflow(
    ppc: dict[str, Any],
    ibr_gens: np.ndarray,
    i_limit: float | None = None,
    tol: float = 1e-8,
    max_outer: int = 300,
    relax: float = 0.7,
) -> IBRPowerFlowResult:
    """Power flow in which IBR reactive output respects the converter current limit.

    A voltage-controlled IBR plant holds its scheduled voltage exactly like a PV
    bus -- until the converter current limit binds. So the algorithm is PV/PQ
    switching, with one crucial difference from the classical version: the limit
    it switches on is not a constant from the case file, it is

        Qmax(V, P) = sqrt( (V * Ilim * Srating)^2 - P^2 )

    which moves every iteration because it depends on the solved voltage. That
    single substitution is the whole experiment: everything else is held fixed,
    so any divergence between this and `classical_powerflow` is attributable to
    the current-limit model and nothing else.

    Outer loop:
      1. Solve with the present PV/PQ assignment.
      2. At each IBR bus compute Q and the capability Qmax(V, P).
      3. Buses exceeding capability are pinned to PQ at +/-Qmax; buses that come
         back inside capability are released to PV again.
      4. Repeat until the binding set and the reactive schedule both stop moving.

    Releasing buses (step 3) matters -- without it the loop can latch a bus at a
    limit it no longer violates, which understates voltage and overstates the error.
    """
    cfg = scenarios()["ibr"]
    # steady state uses the CONTINUOUS rating; 1.2 pu is a transient limit that
    # a converter may only hold for a few cycles, so using it here would credit
    # the plant with reactive support it cannot sustain
    if i_limit is None:
        i_limit = cfg.get("current_limit_continuous_pu", cfg["current_limit_pu"])
    mode = cfg["volt_control"].get("mode", "voltage_regulating")

    work = {**ppc, "bus": ppc["bus"].copy(), "gen": ppc["gen"].copy()}
    ibr_gens = np.asarray(ibr_gens, dtype=int)
    idx = bus_index(work)
    base = work["baseMVA"]

    ibr_bus_rows = np.array(
        sorted({idx[int(work["gen"][g, GEN_BUS])] for g in ibr_gens}), dtype=int,
    )
    ibr_bus_rows = np.array(
        [i for i in ibr_bus_rows if work["bus"][i, 1] != REF], dtype=int,
    )

    # per-IBR-bus active power, MVA rating and scheduled voltage, per-unit
    p_pu, s_rating_pu, v_sched, gen_rows = {}, {}, {}, {}
    for i in ibr_bus_rows:
        key = int(i)
        rows = [g for g in ibr_gens if idx[int(work["gen"][g, GEN_BUS])] == i]
        gen_rows[key] = rows
        p_pu[key] = float(work["gen"][rows, PG].sum()) / base
        s_rating_pu[key] = inverter_rating_pu(ppc, rows)
        # inherit the displaced machine's voltage schedule: an IBR replacing a
        # unit is commissioned to regulate the same point, not to 1.0 pu
        vg = float(ppc["gen"][rows[0], VG])
        v_sched[key] = vg if np.isfinite(vg) and vg > 0 else 1.0

    q_classical = _classical_q_at_buses(ppc, ibr_bus_rows)

    # The ONLY difference from the classical solve is the value written into the
    # IBR generators' QMAX/QMIN. Everything else -- PV/PQ switching, Q-limit
    # enforcement on the remaining synchronous machines -- is handled by the same
    # Newton-Raphson code path. That is what keeps this a controlled experiment:
    # at 0% penetration the two are the same computation, bit for bit.
    v = None
    inner_total = 0
    converged = False
    it = 0
    q_cap_prev: dict[int, float] = {}
    pf = None

    while it < max_outer:
        it += 1

        q_cap: dict[int, float] = {}
        for i in ibr_bus_rows:
            key = int(i)
            vm = float(abs(v[i])) if v is not None else float(v_sched[key])
            cap = q_capability_pu(vm, p_pu[key], i_limit, s_rating_pu[key])
            q_cap[key] = cap
            rows = gen_rows[key]
            work["gen"][rows, QMAX] = 0.0
            work["gen"][rows, QMIN] = 0.0
            work["gen"][rows[0], QMAX] = cap * base
            work["gen"][rows[0], QMIN] = -cap * base
            work["gen"][rows, VG] = v_sched[key]

        pf = solvers.newton_raphson(
            work, tol=1e-11, max_iter=50, enforce_q_limits=True, v0=v,
        )
        if not pf.converged and v is not None:
            pf = solvers.newton_raphson(
                work, tol=1e-11, max_iter=50, enforce_q_limits=True,
            )
        inner_total += pf.iterations
        if not pf.converged:
            return IBRPowerFlowResult(
                False, it, inner_total, pf.v if v is None else v,
                penetration_pct=actual_penetration_pct(ppc, ibr_gens),
                note="inner NR diverged from both warm and flat start",
            )
        v = pf.v

        # capability depends on the solved voltage, so iterate until it settles
        moved = max(
            (abs(q_cap[k] - q_cap_prev.get(k, float("inf"))) for k in q_cap),
            default=0.0,
        )
        q_cap_prev = q_cap
        if moved < tol:
            converged = True
            break

    ybus = yb.build(work)
    s_bus = v * np.conj(ybus @ v)
    q_report = {}
    for i in ibr_bus_rows:
        q_gen = float(s_bus[i].imag) + float(work["bus"][i, 3]) / base   # + QD
        q_report[int(ppc["bus"][i, 0])] = q_gen * base

    # Which IBRs actually hit their converter limit. Taken from the solver's own
    # PV->PQ switching record rather than re-derived from the solution, so it
    # cannot disagree with what the solve actually did.
    ibr_bus_numbers = {int(ppc["bus"][i, 0]) for i in ibr_bus_rows}
    limited = sorted(ibr_bus_numbers.intersection(pf.q_limits_hit))

    return IBRPowerFlowResult(
        converged=converged,
        outer_iterations=it,
        inner_iterations=inner_total,
        v=v if v is not None else np.ones(ppc["bus"].shape[0], dtype=complex),
        q_ibr_mvar=q_report,
        q_classical_mvar=q_classical,
        limited_buses=sorted(limited),
        penetration_pct=actual_penetration_pct(ppc, ibr_gens),
    )


def _classical_q_at_buses(ppc: dict[str, Any], bus_rows: np.ndarray) -> dict[int, float]:
    """Reactive output the classical PV-bus model assigns at the given buses."""
    pf = solvers.newton_raphson(ppc, tol=1e-11)
    if not pf.converged:
        return {}
    ybus = yb.build(ppc)
    s = pf.v * np.conj(ybus @ pf.v)
    out = {}
    for i in bus_rows:
        q_net = float(s[i].imag) * ppc["baseMVA"]
        out[int(ppc["bus"][i, 0])] = q_net + float(ppc["bus"][i, 3])   # + QD
    return out


def classical_powerflow(ppc: dict[str, Any]) -> solvers.PFResult:
    """Baseline: unmodified classical load flow, PV buses with fixed Q limits."""
    return solvers.newton_raphson(ppc, tol=1e-11, enforce_q_limits=True)
