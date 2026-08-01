"""MATPOWER-format case access and column constants.

The "ppc" (MATPOWER case struct) is this project's tool-neutral internal
representation. pandapower produces one; MATPOWER consumes one; our own
solvers, fault code and interchange exporters all read it. Working in ppc
means MATLAB and Python are solving *literally the same numbers*, which is
what makes the cross-tool comparison in WP1 honest rather than approximate.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np

# --- MATPOWER bus matrix columns --------------------------------------------
BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN = range(13)

# --- MATPOWER gen matrix columns --------------------------------------------
GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN = range(10)

# --- MATPOWER branch matrix columns -----------------------------------------
(F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C,
 TAP, SHIFT, BR_STATUS, ANGMIN, ANGMAX) = range(13)

# --- bus types ---------------------------------------------------------------
PQ, PV, REF, NONE = 1, 2, 3, 4

BUS_TYPE_NAME = {PQ: "PQ", PV: "PV", REF: "REF (slack)", NONE: "isolated"}


def load_ppc(case: str) -> dict[str, Any]:
    """Build a MATPOWER-format case dict for one of the IEEE test systems.

    Sourced from pandapower.networks so Python needs no MATLAB to run, but the
    resulting arrays are the standard MATPOWER matrices. A power flow is run
    once to materialise the internal ppc, then bus voltages are reset to flat
    start so downstream solvers are not handed the answer.
    """
    import pandapower as pp
    import pandapower.networks as pn

    if not hasattr(pn, case):
        raise KeyError(f"pandapower.networks has no case {case!r}")

    net = getattr(pn, case)()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pp.runpp(net, numba=False)

    src = net._ppc
    n_line, n_trafo = len(net.line), len(net.trafo)
    n_branch = src["branch"].shape[0]

    # pandapower orders ppc branches as [lines..., transformers...]. Recording
    # the mask here beats guessing from the TAP column, which pandapower fills
    # with 1.0 for lines rather than MATPOWER's 0.0.
    is_xfmr = np.zeros(n_branch, dtype=bool)
    if n_trafo and n_line + n_trafo <= n_branch:
        is_xfmr[n_line: n_line + n_trafo] = True

    ppc: dict[str, Any] = {
        "baseMVA": float(src["baseMVA"]),
        "bus": np.array(src["bus"], dtype=float),
        "gen": np.array(src["gen"], dtype=float),
        "branch": np.array(src["branch"], dtype=float),
        "case": case,
        "is_transformer": is_xfmr,
        # keep the reference Ybus so ybus.build() can be validated against it
        "_ref_Ybus": src["internal"]["Ybus"].toarray().copy(),
        "_net": net,
    }
    _renumber_buses_one_based(ppc)
    ppc["is_slack_gen"] = _slack_gen_mask(ppc)
    ppc["gen_mva_base"] = _gen_mva_base(ppc)
    # flat start: our solvers must find the solution, not be given it
    ppc["bus"][:, VM] = 1.0
    ppc["bus"][:, VA] = 0.0
    return ppc


def _renumber_buses_one_based(ppc: dict[str, Any]) -> None:
    """Renumber buses to 1..n in place, matching IEEE / MATPOWER convention.

    pandapower emits 0-based bus IDs, but every IEEE test-system reference (and
    every fault location in config/scenarios.yaml) uses 1-based numbering. Doing
    the mapping once here means bus 14 of the IEEE 14-bus system is bus 14
    everywhere in this project -- including in the MATLAB deliverables, which
    read the MATPOWER cases natively and are 1-based already.
    """
    old = ppc["bus"][:, BUS_I].astype(int)
    mapping = {int(o): i + 1 for i, o in enumerate(sorted(old.tolist()))}

    ppc["bus"][:, BUS_I] = [mapping[int(b)] for b in old]
    ppc["gen"][:, GEN_BUS] = [mapping[int(b)] for b in ppc["gen"][:, GEN_BUS]]
    ppc["branch"][:, F_BUS] = [mapping[int(b)] for b in ppc["branch"][:, F_BUS]]
    ppc["branch"][:, T_BUS] = [mapping[int(b)] for b in ppc["branch"][:, T_BUS]]
    ppc["bus_renumbering"] = mapping


def _slack_gen_mask(ppc: dict[str, Any]) -> np.ndarray:
    """Which generator rows sit at a REF (slack) bus."""
    ref_buses = set(ppc["bus"][ppc["bus"][:, BUS_TYPE] == REF, BUS_I].astype(int).tolist())
    return np.array([int(b) in ref_buses for b in ppc["gen"][:, GEN_BUS]], dtype=bool)


def _gen_mva_base(ppc: dict[str, Any]) -> np.ndarray:
    """Per-generator MVA base, inferred because the IEEE cases do not carry one.

    pandapower leaves MBASE as NaN, and the slack generator has PMAX = 1e9
    (an unbounded external grid), so neither column can be used directly. The
    rule is declared in config/scenarios.yaml and applied here; slack rows get
    NaN because they are handled separately, via an assumed short-circuit level.
    """
    from .config import scenarios

    md = scenarios()["machine_data"]
    base = ppc["baseMVA"]
    gen = ppc["gen"]
    slack = ppc["is_slack_gen"]

    mbase = gen[:, MBASE].astype(float).copy()
    bad = ~np.isfinite(mbase) | (mbase <= 1.0)      # pandapower writes NaN or 1.0

    mode = md.get("mva_base_from", "pmax")
    if mode == "pmax":
        cand = gen[:, PMAX].astype(float)
        cand = np.where(np.isfinite(cand) & (cand > 0) & (cand < 1e6), cand, np.nan)
    elif mode == "pg":
        cand = np.abs(gen[:, PG].astype(float)) * float(md.get("pg_to_mva_factor", 1.2))
    else:
        cand = np.full(gen.shape[0], base)

    mbase = np.where(bad, cand, mbase)
    # last resort for a finite non-slack machine with no usable rating
    mbase = np.where(np.isfinite(mbase), mbase, np.abs(gen[:, PG]) * 1.2)
    mbase = np.where(np.isfinite(mbase) & (mbase > 0), mbase, base)
    mbase[slack] = np.nan                            # handled via slack SC level
    return mbase


def slack_thevenin_pu(ppc: dict[str, Any]) -> complex:
    """Thevenin impedance of the external grid at the slack bus, per-unit.

    Derived from an assumed short-circuit level rather than a machine reactance,
    because the slack in these cases represents the rest of the interconnection,
    not a specific unit. Z = V^2 / Ssc, split by the configured X/R.
    """
    from .config import scenarios

    cfg = scenarios()["machine_data"]["slack"]
    sc_mva = cfg.get("sc_mva")
    if sc_mva in (None, 0):
        sc_mva = float(cfg.get("sc_multiple_of_load", 10.0)) * total_load_mw(ppc)
    sc_pu = float(sc_mva) / ppc["baseMVA"]
    if sc_pu <= 0:
        raise ValueError("slack short-circuit level resolved to zero")

    z_mag = 1.0 / sc_pu
    xr = float(cfg.get("x_over_r", 10.0))
    x = z_mag * xr / np.sqrt(1.0 + xr ** 2)
    r = x / xr
    return complex(r, x)


def bus_index(ppc: dict[str, Any]) -> dict[int, int]:
    """Map external bus number -> internal row index."""
    return {int(b): i for i, b in enumerate(ppc["bus"][:, BUS_I])}


def bus_types(ppc: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (ref, pv, pq) internal index arrays."""
    types = ppc["bus"][:, BUS_TYPE].astype(int)
    return (
        np.flatnonzero(types == REF),
        np.flatnonzero(types == PV),
        np.flatnonzero(types == PQ),
    )


def bus_injections(ppc: dict[str, Any]) -> np.ndarray:
    """Complex scheduled bus injection Sbus in per-unit (generation - load)."""
    base = ppc["baseMVA"]
    n = ppc["bus"].shape[0]
    idx = bus_index(ppc)

    sbus = -(ppc["bus"][:, PD] + 1j * ppc["bus"][:, QD]) / base

    gen = ppc["gen"]
    live = gen[:, GEN_STATUS] > 0
    for row in gen[live]:
        i = idx[int(row[GEN_BUS])]
        sbus[i] += (row[PG] + 1j * row[QG]) / base
    assert sbus.shape == (n,)
    return sbus


def gen_at_bus(ppc: dict[str, Any]) -> dict[int, list[int]]:
    """Internal bus index -> list of in-service generator row indices."""
    idx = bus_index(ppc)
    out: dict[int, list[int]] = {}
    for g, row in enumerate(ppc["gen"]):
        if row[GEN_STATUS] <= 0:
            continue
        out.setdefault(idx[int(row[GEN_BUS])], []).append(g)
    return out


def total_load_mw(ppc: dict[str, Any]) -> float:
    return float(ppc["bus"][:, PD].sum())


def total_gen_mw(ppc: dict[str, Any]) -> float:
    gen = ppc["gen"]
    return float(gen[gen[:, GEN_STATUS] > 0, PG].sum())


def summary(ppc: dict[str, Any]) -> str:
    ref, pv, pq = bus_types(ppc)
    return (
        f"{ppc['case']}: {ppc['bus'].shape[0]} bus "
        f"({len(ref)} slack / {len(pv)} PV / {len(pq)} PQ), "
        f"{int((ppc['gen'][:, GEN_STATUS] > 0).sum())} gen, "
        f"{ppc['branch'].shape[0]} branch, "
        f"load {total_load_mw(ppc):.1f} MW"
    )
