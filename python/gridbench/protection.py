"""Protection elements: mho distance, negative-sequence directional, and inverse-time
overcurrent -- evaluated against classical and IBR-aware fault solutions.

This is where the WP3 numbers acquire consequences. A 26% drop in fault current
and a 20-degree shift in fault current angle are abstractions until you ask what
a relay does with them:

  * an inverse-time overcurrent element may simply not pick up;
  * a mho distance element measures an apparent impedance that no longer
    corresponds to the distance to the fault, so Zone 1 under- or over-reaches;
  * a negative-sequence directional element has nothing to polarise on when the
    inverter injects no negative-sequence current -- which is precisely why
    IEEE Std 2800-2022 mandates that it does.

Every element here is deliberately textbook. The point of the study is not a
novel relay; it is that STANDARD, correctly-set relays misoperate when the
source behind them stops being a synchronous machine.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np

from .config import scenarios
from .faults import FaultResult, SequenceNetworks, seq_to_phase
from .metrics import wrap_deg
from .ppc import BR_R, BR_X, F_BUS, T_BUS, bus_index


# =============================================================================
# measured quantities at a relay location
# =============================================================================
@dataclass
class RelayMeasurement:
    """Voltages and currents a relay at `from_bus` sees for a fault at `fault_bus`."""

    from_bus: int
    to_bus: int
    fault_bus: int
    v_phase: np.ndarray          # 3-vector, per-unit
    i_phase: np.ndarray          # 3-vector, per-unit
    v1: complex
    v2: complex
    v0: complex
    i1: complex
    i2: complex
    i0: complex

    @property
    def i2_magnitude(self) -> float:
        return float(abs(self.i2))


def measure_at_relay(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    result: FaultResult,
    branch: int,
) -> RelayMeasurement:
    """Compute what a relay at the FROM end of `branch` measures during the fault.

    Line current is obtained from the sequence voltage profiles across the
    branch series impedance, one sequence at a time, then transformed to phase
    quantities. Shunt charging is neglected at the relay point, which is
    standard for distance protection.
    """
    br = ppc["branch"]
    idx = bus_index(ppc)
    f = idx[int(br[branch, F_BUS])]
    t = idx[int(br[branch, T_BUS])]

    z1_line = complex(br[branch, BR_R], br[branch, BR_X])
    sd = scenarios()["sequence_data"]
    z0_line = complex(
        br[branch, BR_R] * sd["line_r0_over_r1"],
        br[branch, BR_X] * sd["line_x0_over_x1"],
    )

    i1 = (result.v1[f] - result.v1[t]) / z1_line if abs(z1_line) > 0 else 0j
    i2 = (result.v2[f] - result.v2[t]) / z1_line if abs(z1_line) > 0 else 0j
    i0 = (result.v0[f] - result.v0[t]) / z0_line if abs(z0_line) > 0 else 0j

    return RelayMeasurement(
        from_bus=int(br[branch, F_BUS]),
        to_bus=int(br[branch, T_BUS]),
        fault_bus=result.bus,
        v_phase=seq_to_phase(result.v1[f], result.v2[f], result.v0[f]),
        i_phase=seq_to_phase(i1, i2, i0),
        v1=complex(result.v1[f]), v2=complex(result.v2[f]), v0=complex(result.v0[f]),
        i1=complex(i1), i2=complex(i2), i0=complex(i0),
    )


# =============================================================================
# mho distance element
# =============================================================================
@dataclass
class DistanceVerdict:
    apparent_z: complex
    true_z: complex
    zone: int                    # 0 = no pickup, 1/2/3 = tripping zone
    operate_time_s: float
    reach_error_pct: float

    @property
    def tripped(self) -> bool:
        return self.zone > 0


def apparent_impedance(m: RelayMeasurement, fault_type: str) -> complex:
    """Apparent impedance seen by the correct measuring loop for the fault type.

    Phase loops use (Va - Vb) / (Ia - Ib); ground loops use Va / (Ia + k0*3I0)
    with the residual compensation factor k0 = (Z0 - Z1) / (3 Z1). Selecting the
    right loop matters: this is exactly the "phase selection" function that IBR
    fault current is known to confuse.
    """
    v, i = m.v_phase, m.i_phase
    if fault_type in ("LL", "LLG"):
        # b-c loop
        denom = i[1] - i[2]
        return (v[1] - v[2]) / denom if abs(denom) > 1e-12 else complex(np.inf, np.inf)
    if fault_type == "3LG":
        denom = i[0] - i[1]
        return (v[0] - v[1]) / denom if abs(denom) > 1e-12 else complex(np.inf, np.inf)

    # SLG on phase a, with residual compensation
    sd = scenarios()["sequence_data"]
    k0 = (complex(sd["line_r0_over_r1"], 0) - 1.0) / 3.0   # approx (Z0-Z1)/(3Z1)
    denom = i[0] + k0 * 3.0 * m.i0
    return v[0] / denom if abs(denom) > 1e-12 else complex(np.inf, np.inf)


def mho_distance(
    m: RelayMeasurement,
    z_line: complex,
    fault_type: str,
    z_true: complex | None = None,
) -> DistanceVerdict:
    """Three-zone mho distance element.

    A mho characteristic is a circle through the origin whose diameter lies along
    the line angle. The element operates when the apparent impedance falls inside
    the circle, which for a self-polarised mho reduces to

        |Z_apparent - Zr/2| <= |Zr/2|      where Zr is the zone reach.
    """
    cfg = scenarios()["protection"]["distance_relay"]
    z_app = apparent_impedance(m, fault_type)
    z_true = z_line if z_true is None else z_true

    zones = [
        (1, cfg["zone1_reach_pct"], cfg["zone1_delay_s"]),
        (2, cfg["zone2_reach_pct"], cfg["zone2_delay_s"]),
        (3, cfg["zone3_reach_pct"], cfg["zone3_delay_s"]),
    ]

    picked, t_op = 0, float("inf")
    if np.isfinite(z_app.real) and np.isfinite(z_app.imag):
        for zone, reach_pct, delay in zones:
            zr = z_line * (reach_pct / 100.0)
            if abs(z_app - zr / 2.0) <= abs(zr / 2.0) + 1e-12:
                picked, t_op = zone, delay
                break

    reach_err = (
        abs(z_app - z_true) / abs(z_true) * 100.0
        if abs(z_true) > 1e-12 and np.isfinite(abs(z_app)) else float("nan")
    )
    return DistanceVerdict(z_app, z_true, picked, t_op, float(reach_err))


# =============================================================================
# negative-sequence directional element
# =============================================================================
@dataclass
class DirectionalVerdict:
    declared: Literal["forward", "reverse", "no_decision"]
    torque_angle_deg: float
    i2_magnitude: float
    polarisable: bool


def negative_sequence_directional(
    m: RelayMeasurement,
    z1_line: complex,
    i2_threshold: float = 0.05,
) -> DirectionalVerdict:
    """32Q negative-sequence directional element.

    Compares the angle of I2 against -V2 rotated by the line angle. A forward
    fault gives a torque angle near zero.

    `i2_threshold` is the supervising pickup: below it the element refuses to
    make a decision, because polarising on noise is worse than not operating.
    This threshold is the crux of the IEEE 2800 question -- a pre-2800 inverter
    injecting no negative-sequence current leaves the element BLIND, so it
    returns 'no_decision' rather than a wrong answer.
    """
    i2_mag = float(abs(m.i2))
    if i2_mag < i2_threshold or abs(m.v2) < 1e-6:
        return DirectionalVerdict("no_decision", float("nan"), i2_mag, False)

    line_angle = np.angle(z1_line)
    torque = wrap_deg(np.degrees(np.angle(m.i2) - np.angle(-m.v2) + line_angle))
    declared = "forward" if abs(torque) <= 90.0 else "reverse"
    return DirectionalVerdict(declared, float(torque), i2_mag, True)


# =============================================================================
# inverse-time overcurrent
# =============================================================================
def iec_operating_time(
    i_measured: float, i_pickup: float, tms: float = 0.1, curve: str = "IEC_standard_inverse",
) -> float:
    """IEC 60255 inverse-time characteristic. Returns inf if the element does not pick up.

    An inverter-fed fault frequently fails to reach pickup at all, which is the
    simplest and most consequential protection failure in the whole study: not a
    mis-timed trip, but no trip.
    """
    constants = {
        "IEC_standard_inverse": (0.14, 0.02),
        "IEC_very_inverse": (13.5, 1.0),
        "IEC_extremely_inverse": (80.0, 2.0),
    }
    k, alpha = constants.get(curve, constants["IEC_standard_inverse"])
    if i_pickup <= 0:
        return float("inf")
    m = i_measured / i_pickup
    if m <= 1.0:
        return float("inf")            # no pickup
    return float(tms * k / (m ** alpha - 1.0))


# =============================================================================
# study driver
# =============================================================================
@dataclass
class ProtectionOutcome:
    branch: int
    from_bus: int
    to_bus: int
    fault_bus: int
    fault_type: str
    method: str
    zone: int
    operate_time_s: float
    reach_error_pct: float
    directional: str
    torque_angle_deg: float
    i2_magnitude: float
    oc_operate_time_s: float
    picked_up: bool


def evaluate(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    result: FaultResult,
    branch: int,
    i_pickup_pu: float,
) -> ProtectionOutcome:
    """Run all three elements at one relay location for one fault solution."""
    br = ppc["branch"]
    z_line = complex(br[branch, BR_R], br[branch, BR_X])
    m = measure_at_relay(ppc, seq, result, branch)

    dist = mho_distance(m, z_line, result.fault_type)
    direc = negative_sequence_directional(m, z_line)
    i_max = float(np.max(np.abs(m.i_phase)))
    t_oc = iec_operating_time(i_max, i_pickup_pu)

    return ProtectionOutcome(
        branch=branch, from_bus=m.from_bus, to_bus=m.to_bus,
        fault_bus=result.bus, fault_type=result.fault_type, method=result.method,
        zone=dist.zone, operate_time_s=dist.operate_time_s,
        reach_error_pct=dist.reach_error_pct,
        directional=direc.declared, torque_angle_deg=direc.torque_angle_deg,
        i2_magnitude=direc.i2_magnitude,
        oc_operate_time_s=t_oc, picked_up=np.isfinite(t_oc),
    )


def misoperation_flags(
    classical: ProtectionOutcome, ibr: ProtectionOutcome,
) -> dict[str, bool]:
    """Compare an IBR-aware outcome against the classical one it was set from.

    A relay is commissioned using a classical short-circuit study. These flags
    ask, for each element, whether the settings derived that way still behave
    correctly once the sources behind them are inverters.
    """
    return {
        "zone_changed": classical.zone != ibr.zone,
        "lost_distance_trip": classical.zone > 0 and ibr.zone == 0,
        "zone1_reach_lost": classical.zone == 1 and ibr.zone != 1,
        "directional_lost": classical.directional != "no_decision"
                            and ibr.directional == "no_decision",
        "directional_reversed": (
            classical.directional in ("forward", "reverse")
            and ibr.directional in ("forward", "reverse")
            and classical.directional != ibr.directional
        ),
        "lost_oc_pickup": classical.picked_up and not ibr.picked_up,
        "slower_oc": (
            np.isfinite(classical.oc_operate_time_s)
            and np.isfinite(ibr.oc_operate_time_s)
            and ibr.oc_operate_time_s > classical.oc_operate_time_s * 1.2
        ),
    }
