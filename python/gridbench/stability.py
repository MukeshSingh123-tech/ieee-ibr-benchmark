"""Voltage stability and network sensitivity: continuation power flow, L-index,
PTDF/LODF, and N-1 contingency screening.

These are the tools that answer "how much margin is left", as opposed to "what
is the operating point". They matter more, not less, with high IBR penetration:
displacing synchronous machines removes reactive reserve, which is what voltage
stability margin is made of, so the same network carries less margin at the same
dispatch. The continuation power flow here quantifies exactly that.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from . import solvers, ybus as yb
from .ppc import (
    BR_STATUS, BR_X, F_BUS, GEN_BUS, GEN_STATUS, PD, PG, QD, T_BUS,
    bus_index, bus_types,
)


# =============================================================================
# continuation power flow
# =============================================================================
@dataclass
class CPFResult:
    lambdas: list[float] = field(default_factory=list)
    vm_traces: list[np.ndarray] = field(default_factory=list)
    lambda_max: float = 0.0
    critical_bus: int = -1
    converged: bool = False

    @property
    def loading_margin_pct(self) -> float:
        """Percentage load increase available before voltage collapse."""
        return (self.lambda_max - 1.0) * 100.0

    def nose_curve(self, bus_row: int) -> tuple[np.ndarray, np.ndarray]:
        """(lambda, |V|) trace for one bus -- the classic PV nose curve."""
        return (np.array(self.lambdas),
                np.array([v[bus_row] for v in self.vm_traces]))


def continuation_power_flow(
    ppc: dict[str, Any],
    step: float = 0.05,
    min_step: float = 1e-4,
    max_points: int = 400,
    enforce_q_limits: bool = True,
) -> CPFResult:
    """Trace the PV nose curve by predictor-corrector continuation.

    Loads (and generation, to keep the balance) are scaled by a continuation
    parameter lambda. Near the nose the Jacobian becomes singular and plain
    Newton-Raphson fails, so the step is halved on failure and retried; the
    largest lambda that still converges is the loading margin.

    This is a simple arc-length-free continuation: adequate for finding the nose
    and cheap to implement, but it cannot trace the lower (unstable) branch. The
    report states that limitation rather than implying a full CPF.
    """
    result = CPFResult()
    base_bus = ppc["bus"].copy()
    base_gen = ppc["gen"].copy()

    lam = 1.0
    v_prev = None
    current_step = step

    while len(result.lambdas) < max_points and current_step >= min_step:
        trial = {**ppc, "bus": base_bus.copy(), "gen": base_gen.copy()}
        trial["bus"][:, PD] = base_bus[:, PD] * lam
        trial["bus"][:, QD] = base_bus[:, QD] * lam
        # scale non-slack generation with load so the slack does not absorb it all
        non_slack = (~ppc["is_slack_gen"]) & (base_gen[:, GEN_STATUS] > 0)
        trial["gen"][non_slack, PG] = base_gen[non_slack, PG] * lam

        pf = solvers.newton_raphson(
            trial, tol=1e-10, max_iter=40,
            enforce_q_limits=enforce_q_limits, v0=v_prev,
        )

        if pf.converged and np.all(np.abs(pf.v) > 0.3):
            result.lambdas.append(lam)
            result.vm_traces.append(pf.vm.copy())
            result.lambda_max = lam
            v_prev = pf.v
            lam += current_step
        else:
            # past the nose (or nearly): refine and try again
            lam -= current_step
            current_step /= 2.0
            lam += current_step

    result.converged = len(result.lambdas) > 1
    if result.vm_traces:
        result.critical_bus = int(ppc["bus"][int(np.argmin(result.vm_traces[-1])), 0])
    return result


# =============================================================================
# L-index voltage stability indicator
# =============================================================================
def l_index(ppc: dict[str, Any], v: np.ndarray) -> np.ndarray:
    """Kessel-Glavitsch L-index per load bus. 0 = no load, 1 = voltage collapse.

        L_j = | 1 - sum_{i in gen} F_ji * V_i / V_j |

    where F = -[Y_LL]^-1 [Y_LG] is derived from the partitioned admittance
    matrix. Unlike the CPF it needs no continuation -- one solved power flow
    gives a proximity-to-collapse figure for every load bus, which is why it is
    used for online monitoring.
    """
    ref, pv, pq = bus_types(ppc)
    gen_rows = np.sort(np.concatenate([ref, pv]))
    load_rows = pq

    n = ppc["bus"].shape[0]
    out = np.zeros(n)
    if load_rows.size == 0 or gen_rows.size == 0:
        return out

    ybus = yb.build(ppc)
    y_ll = ybus[np.ix_(load_rows, load_rows)]
    y_lg = ybus[np.ix_(load_rows, gen_rows)]

    try:
        f_lg = -np.linalg.solve(y_ll, y_lg)
    except np.linalg.LinAlgError:
        out[load_rows] = np.nan
        return out

    for k, j in enumerate(load_rows):
        if abs(v[j]) < 1e-9:
            out[j] = 1.0
            continue
        out[j] = float(abs(1.0 - np.sum(f_lg[k, :] * v[gen_rows]) / v[j]))
    return out


# =============================================================================
# linear sensitivities
# =============================================================================
def ptdf(ppc: dict[str, Any], slack_bus: int | None = None) -> np.ndarray:
    """Power Transfer Distribution Factors, shape (n_branch, n_bus).

    PTDF[l, i] is the change in MW flow on branch l per MW injected at bus i and
    withdrawn at the slack. Derived from the DC power flow: linear, lossless,
    and the workhorse of security-constrained dispatch and market clearing.
    """
    bus, br = ppc["bus"], ppc["branch"]
    n, m = bus.shape[0], br.shape[0]
    idx = bus_index(ppc)

    ref, _, _ = bus_types(ppc)
    slack = int(ref[0]) if slack_bus is None else idx[int(slack_bus)]

    b_series = np.zeros(m)
    f_idx = np.zeros(m, dtype=int)
    t_idx = np.zeros(m, dtype=int)
    for k in range(m):
        f_idx[k] = idx[int(br[k, F_BUS])]
        t_idx[k] = idx[int(br[k, T_BUS])]
        x = br[k, BR_X]
        b_series[k] = 1.0 / x if (br[k, BR_STATUS] > 0 and abs(x) > 1e-12) else 0.0

    # branch-bus incidence and susceptance matrices
    incidence = np.zeros((m, n))
    incidence[np.arange(m), f_idx] = 1.0
    incidence[np.arange(m), t_idx] = -1.0

    bf = np.diag(b_series) @ incidence           # branch flow matrix
    b_bus = incidence.T @ bf                     # nodal susceptance matrix

    non_slack = np.array([i for i in range(n) if i != slack])
    b_red = b_bus[np.ix_(non_slack, non_slack)]

    out = np.zeros((m, n))
    try:
        theta_sens = np.linalg.solve(b_red, np.eye(non_slack.size))
    except np.linalg.LinAlgError:
        return out
    out[:, non_slack] = bf[:, non_slack] @ theta_sens
    return out


def lodf(ppc: dict[str, Any], ptdf_matrix: np.ndarray | None = None) -> np.ndarray:
    """Line Outage Distribution Factors, shape (n_branch, n_branch).

    LODF[l, k] is the fraction of branch k's pre-outage flow that appears on
    branch l after k trips. This gives an entire N-1 screen from one matrix
    instead of one power flow per contingency -- the standard way real-time
    contingency analysis keeps up with the clock.

    The diagonal is set to -1 (an outaged line loses all its own flow), and a
    column is zeroed where the outage would island the network (PTDF self-term
    equal to 1), which the linear theory cannot represent.
    """
    br = ppc["branch"]
    m = br.shape[0]
    idx = bus_index(ppc)
    h = ptdf(ppc) if ptdf_matrix is None else ptdf_matrix

    out = np.zeros((m, m))
    for k in range(m):
        f = idx[int(br[k, F_BUS])]
        t = idx[int(br[k, T_BUS])]
        self_term = h[k, f] - h[k, t]
        if abs(1.0 - self_term) < 1e-9:          # outage islands the system
            out[:, k] = 0.0
            out[k, k] = -1.0
            continue
        out[:, k] = (h[:, f] - h[:, t]) / (1.0 - self_term)
        out[k, k] = -1.0
    return out


# =============================================================================
# N-1 contingency screening
# =============================================================================
@dataclass
class ContingencyResult:
    branch: int
    from_bus: int
    to_bus: int
    converged: bool
    max_loading_pct: float
    worst_branch: int
    min_vm_pu: float
    max_vm_pu: float
    n_violations: int
    islanded: bool = False

    @property
    def status(self) -> str:
        if self.islanded:
            return "islands network"
        return "solved" if self.converged else "diverged"


def islands_network(ppc: dict[str, Any], branch: int) -> bool:
    """Would removing `branch` disconnect the network?

    Worth separating from divergence: an outage that islands buses is a
    topological fact the operator already knows about, whereas a converged-but-
    violating case or a genuine solver divergence are different findings. Lumping
    them together as "did not converge" would overstate how many contingencies
    are numerically troublesome.
    """
    br = ppc["branch"]
    idx = bus_index(ppc)
    n = ppc["bus"].shape[0]

    adjacency: list[set[int]] = [set() for _ in range(n)]
    for k in range(br.shape[0]):
        if k == branch or br[k, BR_STATUS] <= 0:
            continue
        f, t = idx[int(br[k, F_BUS])], idx[int(br[k, T_BUS])]
        adjacency[f].add(t)
        adjacency[t].add(f)

    seen = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for nb in adjacency[node]:
            if nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return len(seen) < n


def _violation_set(
    ppc: dict[str, Any], pf, vm_limits: tuple[float, float], rate_mva: float | None,
) -> set[str]:
    """Identify which limits are violated, as a set of labels."""
    bad = set()
    for i, vm in enumerate(pf.vm):
        bus = int(ppc["bus"][i, 0])
        if vm < vm_limits[0]:
            bad.add(f"vlow:{bus}")
        elif vm > vm_limits[1]:
            bad.add(f"vhigh:{bus}")

    sf, st = yb.branch_flows(ppc, pf.v)
    flow_mva = np.maximum(np.abs(sf), np.abs(st)) * ppc["baseMVA"]
    rating = ppc["branch"][:, 5].copy()
    rating[rating <= 0] = rate_mva if rate_mva else np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        loading = flow_mva / rating * 100.0
    for k, load in enumerate(loading):
        if np.isfinite(load) and load > 100.0:
            bad.add(f"thermal:{k}")
    return bad


def n1_screen(
    ppc: dict[str, Any],
    rate_mva: float | None = None,
    vm_limits: tuple[float, float] = (0.95, 1.05),
) -> list[ContingencyResult]:
    """Full AC N-1 branch-outage screen.

    Slower than the LODF screen but exact, including reactive power and voltage
    limits that the DC approximation cannot see.

    Violations are counted RELATIVE TO THE PRE-CONTINGENCY CASE -- that is, only
    limits the outage newly breaks. Counting absolute excursions instead makes
    the ranking useless: several IEEE cases sit slightly outside 0.95-1.05 pu at
    some bus even intact, so every contingency scores as "violating" and nothing
    is distinguishable from anything else. New violations are also what an
    operator actually acts on.
    """
    br = ppc["branch"]
    m = br.shape[0]
    out: list[ContingencyResult] = []

    base_pf = solvers.newton_raphson(ppc, tol=1e-9, max_iter=40, enforce_q_limits=True)
    base_bad = (_violation_set(ppc, base_pf, vm_limits, rate_mva)
                if base_pf.converged else set())

    for k in range(m):
        if br[k, BR_STATUS] <= 0:
            continue

        if islands_network(ppc, k):
            out.append(ContingencyResult(
                k, int(br[k, F_BUS]), int(br[k, T_BUS]), False,
                float("nan"), -1, float("nan"), float("nan"), 0, islanded=True,
            ))
            continue

        trial = {**ppc, "branch": br.copy()}
        trial["branch"][k, BR_STATUS] = 0

        try:
            pf = solvers.newton_raphson(trial, tol=1e-9, max_iter=40,
                                        enforce_q_limits=True)
        except np.linalg.LinAlgError:
            pf = None

        if pf is None or not pf.converged:
            out.append(ContingencyResult(
                k, int(br[k, F_BUS]), int(br[k, T_BUS]), False,
                float("nan"), -1, float("nan"), float("nan"), -1,
            ))
            continue

        sf, st = yb.branch_flows(trial, pf.v)
        flow_mva = np.maximum(np.abs(sf), np.abs(st)) * ppc["baseMVA"]
        rating = br[:, 5].copy()
        rating[rating <= 0] = rate_mva if rate_mva else np.nan
        with np.errstate(invalid="ignore", divide="ignore"):
            loading = flow_mva / rating * 100.0
        loading[k] = np.nan                       # the outaged branch carries nothing

        worst = int(np.nanargmax(loading)) if np.any(np.isfinite(loading)) else -1
        # only limits this outage NEWLY breaks
        n_viol = len(_violation_set(trial, pf, vm_limits, rate_mva) - base_bad)

        out.append(ContingencyResult(
            k, int(br[k, F_BUS]), int(br[k, T_BUS]), True,
            float(np.nanmax(loading)) if worst >= 0 else float("nan"),
            worst, float(pf.vm.min()), float(pf.vm.max()), n_viol,
        ))

    # rank worst-first: genuine divergence, then islanding, then by violations
    def rank(c: ContingencyResult) -> tuple[int, float, float]:
        severity = 0 if (not c.converged and not c.islanded) else (1 if c.islanded else 2)
        return (severity,
                -(c.n_violations if c.converged else 0),
                -(c.max_loading_pct if np.isfinite(c.max_loading_pct) else 0.0))

    out.sort(key=rank)
    return out
