"""System strength metrics: SCR, WSCR, CSCR, ESCR.

System strength is the single most important screening quantity in modern
interconnection studies, and it has no analogue in classical power systems
teaching -- because with synchronous machines it was never scarce.

The metrics differ in what they account for:

  SCR   Short-Circuit Ratio at one point of interconnection. Ignores the fact
        that neighbouring IBRs consume each other's strength.
  WSCR  Weighted SCR (ERCOT). Aggregates a GROUP of IBRs that interact, and is
        the metric that first exposed multi-inverter instability in West Texas.
  CSCR  Composite SCR. Group treated through one equivalent impedance.
  ESCR  Equivalent SCR. Discounts reactive compensation, which inflates fault
        level without contributing to synchronising strength.

The usual interconnection screen is SCR >= 3.0. Below that, grid-following
converters risk PLL synchronisation instability and sub-synchronous control
interaction -- and, critically for this project, phasor-domain tools cannot see
either phenomenon. Low SCR is therefore the flag that says "your load flow
answer is no longer sufficient evidence".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .config import scenarios
from .faults import SequenceNetworks, short_circuit_mva
from .ppc import GEN_BUS, PG, bus_index


@dataclass
class StrengthResult:
    bus: int
    sc_mva: float
    p_ibr_mw: float
    scr: float
    classification: str

    @property
    def is_weak(self) -> bool:
        return self.scr < scenarios()["system_strength"]["weak_grid_threshold"]


def classify(scr: float) -> str:
    """Industry-conventional bands for interpreting SCR."""
    if not np.isfinite(scr):
        return "no IBR at bus"
    if scr >= 5.0:
        return "strong"
    if scr >= 3.0:
        return "moderate"
    if scr >= 2.0:
        return "weak"
    return "very weak"


def scr_at_bus(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    bus: int,
    p_ibr_mw: float,
) -> StrengthResult:
    """Plain SCR = short-circuit MVA / IBR active power rating at one bus.

    Note the short-circuit MVA must be computed with the IBR's own contribution
    EXCLUDED -- an inverter cannot supply its own system strength. That exclusion
    is handled upstream by `faults.build_sequence_networks`, which omits machine
    shunts at IBR buses.
    """
    ssc = short_circuit_mva(ppc, seq, bus)
    scr = float(ssc / p_ibr_mw) if p_ibr_mw > 1e-9 else float("inf")
    return StrengthResult(int(bus), ssc, p_ibr_mw, scr, classify(scr))


def weighted_scr(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    ibr_bus_mw: dict[int, float],
) -> float:
    """ERCOT Weighted Short-Circuit Ratio for an interacting group of IBRs.

                    sum_i ( Ssc_i * P_i )
        WSCR = -----------------------------
                    ( sum_i P_i )^2

    The squared denominator is the whole point: it penalises concentration. Ten
    plants of 100 MW clustered in one weak pocket score far worse than the same
    1000 MW spread across strong buses, which plain per-bus SCR would miss.
    """
    if not ibr_bus_mw:
        return float("inf")
    total_mw = sum(ibr_bus_mw.values())
    if total_mw <= 1e-9:
        return float("inf")
    numerator = sum(
        short_circuit_mva(ppc, seq, bus) * mw for bus, mw in ibr_bus_mw.items()
    )
    return float(numerator / total_mw ** 2)


def composite_scr(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    ibr_bus_mw: dict[int, float],
) -> float:
    """Composite SCR: the group seen through one equivalent Thevenin impedance.

    The equivalent impedance is the MW-weighted mean of the diagonal Thevenin
    impedances, which credits the group for being electrically distributed
    without assuming the buses are independent.
    """
    if not ibr_bus_mw:
        return float("inf")
    total_mw = sum(ibr_bus_mw.values())
    if total_mw <= 1e-9:
        return float("inf")

    idx = bus_index(ppc)
    z_eq = sum(
        seq.z1[idx[int(bus)], idx[int(bus)]] * mw for bus, mw in ibr_bus_mw.items()
    ) / total_mw
    ssc_eq = ppc["baseMVA"] / abs(z_eq) if abs(z_eq) > 1e-12 else float("inf")
    return float(ssc_eq / total_mw)


def equivalent_scr(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    ibr_bus_mw: dict[int, float],
    q_compensation_mvar: float = 0.0,
) -> float:
    """ESCR: SCR discounted by reactive compensation at the point of connection.

    Shunt capacitors and SVCs raise measured fault level without adding
    synchronising torque, so counting them as "strength" overstates the grid.
    ESCR subtracts them -- which is why a plant can pass an SCR screen and still
    fail an ESCR one.
    """
    if not ibr_bus_mw:
        return float("inf")
    total_mw = sum(ibr_bus_mw.values())
    if total_mw <= 1e-9:
        return float("inf")
    ssc = sum(short_circuit_mva(ppc, seq, bus) for bus in ibr_bus_mw)
    return float((ssc - q_compensation_mvar) / total_mw)


def profile(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    ibr_gens: np.ndarray,
) -> dict[str, Any]:
    """Full system-strength picture for one penetration scenario."""
    idx_map = {}
    for g in np.asarray(ibr_gens, dtype=int):
        bus = int(ppc["gen"][g, GEN_BUS])
        idx_map[bus] = idx_map.get(bus, 0.0) + float(ppc["gen"][g, PG])

    per_bus = {
        bus: scr_at_bus(ppc, seq, bus, mw) for bus, mw in idx_map.items()
    }
    scrs = [r.scr for r in per_bus.values() if np.isfinite(r.scr)]

    return {
        "per_bus": per_bus,
        "wscr": weighted_scr(ppc, seq, idx_map),
        "cscr": composite_scr(ppc, seq, idx_map),
        "escr": equivalent_scr(ppc, seq, idx_map),
        "min_scr": float(min(scrs)) if scrs else float("inf"),
        "mean_scr": float(np.mean(scrs)) if scrs else float("inf"),
        "n_weak_buses": sum(1 for r in per_bus.values() if r.is_weak),
        "total_ibr_mw": sum(idx_map.values()),
    }
