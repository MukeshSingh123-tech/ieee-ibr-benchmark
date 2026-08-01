"""Short-circuit analysis: classical symmetrical components, and an IBR-aware solver.

Two solvers live here, and the difference between them is the point of WP3.

`solve_fault_classical`
    Textbook method. Synchronous machines are voltage sources behind X"d, the
    three sequence networks are independent, superposition holds, and the answer
    is a closed-form expression. This is what every commercial short-circuit tool
    (and every power systems course) does.

`solve_fault_ibr_aware`
    Inverter-based resources are NOT voltage sources. They are current-limited,
    controlled sources whose output depends on the very terminal voltage the
    fault produces. Two consequences break the classical method outright:

      1. The network becomes NONLINEAR -- the injection depends on the solution,
         so the problem must be iterated rather than solved in closed form.
      2. Negative-sequence current is a CONTROL CHOICE (IEEE Std 2800-2022
         requires I2 = K2*V2), not a machine property. The sequence networks
         stop being independent, so the classical "open-circuit voltage exists
         only in the positive sequence" assumption fails.

Both solvers share the same boundary-condition formulation below, which is
written in terms of open-circuit sequence voltages (V1_oc, V2_oc, V0_oc) at the
fault bus. The classical case is simply the special case V2_oc = V0_oc = 0.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from . import ybus as yb
from .config import scenarios
from .ppc import (
    BR_B, BR_R, BR_STATUS, BR_X, PD, PG, QD,
    bus_index, gen_at_bus, slack_thevenin_pu,
)

FaultType = Literal["3LG", "SLG", "LL", "LLG"]

# symmetrical-component transformation, A = [[1,1,1],[1,a^2,a],[1,a,a^2]]
A_OP = np.exp(2j * np.pi / 3)
A_MATRIX = np.array([
    [1, 1, 1],
    [1, A_OP ** 2, A_OP],
    [1, A_OP, A_OP ** 2],
], dtype=complex)


# =============================================================================
# sequence network construction
# =============================================================================
@dataclass
class SequenceNetworks:
    """Positive/negative/zero sequence bus admittance and impedance matrices."""

    y1: np.ndarray
    y2: np.ndarray
    y0: np.ndarray
    z1: np.ndarray
    z2: np.ndarray
    z0: np.ndarray
    e_internal: np.ndarray          # voltage-source Norton current injections (pos seq)
    ibr_buses: np.ndarray           # internal indices modelled as IBRs (GFL + GFM)
    sync_buses: np.ndarray          # internal indices with synchronous machines
    v_prefault: np.ndarray          # prefault bus voltages from the load flow
    gfm_buses: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    gfl_buses: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    @property
    def n(self) -> int:
        return self.y1.shape[0]

    @property
    def voltage_source_buses(self) -> np.ndarray:
        """Buses that impose a voltage rather than inject a current.

        Synchronous machines AND grid-forming converters. This is the set that
        determines whether the phasor fault problem is well-posed at all: with
        none of them, every source is a current injection whose magnitude
        depends on a voltage nothing establishes.
        """
        return np.array(sorted(set(self.sync_buses.tolist())
                               | set(self.gfm_buses.tolist())), dtype=int)


def _machine_shunt(ppc: dict[str, Any], x_pu: float, r_over_x: float,
                   buses: np.ndarray, seq: str = "pos") -> np.ndarray:
    """Diagonal shunt admittance from machine reactance, on SYSTEM base.

    Slack generators are excluded from the reactance model and given the
    external-grid Thevenin impedance instead (see `ppc.slack_thevenin_pu`), since
    an unbounded external grid has no meaningful machine rating. In the zero
    sequence the external grid is assumed solidly grounded with the same Z.
    """
    n = ppc["bus"].shape[0]
    ysh = np.zeros(n, dtype=complex)
    gens = gen_at_bus(ppc)
    mva_base = ppc["gen_mva_base"]
    is_slack = ppc["is_slack_gen"]
    z_slack = slack_thevenin_pu(ppc)

    for i in buses:
        for g in gens.get(int(i), []):
            if is_slack[g]:
                ysh[i] += 1.0 / z_slack
                continue
            mbase = mva_base[g]
            if not np.isfinite(mbase) or mbase <= 0:
                mbase = ppc["baseMVA"]
            x_sys = x_pu * ppc["baseMVA"] / mbase      # machine base -> system base
            ysh[i] += 1.0 / complex(r_over_x * x_sys, x_sys)
    return ysh


def _load_shunt(ppc: dict[str, Any], v_pre: np.ndarray) -> np.ndarray:
    """Constant-impedance load representation from the prefault solution."""
    cfg = scenarios()["sequence_data"]
    if cfg.get("load_model") == "neglect":
        return np.zeros(ppc["bus"].shape[0], dtype=complex)
    s_load = (ppc["bus"][:, PD] + 1j * ppc["bus"][:, QD]) / ppc["baseMVA"]
    vm2 = np.abs(v_pre) ** 2
    vm2[vm2 < 1e-9] = 1.0
    return np.conj(s_load) / vm2          # Y = conj(S)/|V|^2


def _zero_sequence_ppc(ppc: dict[str, Any]) -> dict[str, Any]:
    """Derive a zero-sequence network from positive-sequence data.

    IEEE test cases have no zero-sequence data, so it is synthesised with the
    typical multipliers declared in config/scenarios.yaml. Transformer branches
    (TAP != 0) are opened, standing in for the delta-wye(g) connection that
    blocks zero-sequence current -- a coarse but standard approximation.
    """
    cfg = scenarios()["sequence_data"]
    p0 = {**ppc, "bus": ppc["bus"].copy(), "branch": ppc["branch"].copy()}
    p0["branch"][:, BR_R] *= cfg["line_r0_over_r1"]
    p0["branch"][:, BR_X] *= cfg["line_x0_over_x1"]
    p0["branch"][:, BR_B] *= cfg["line_b0_over_b1"]
    if cfg.get("transformer_blocks_zero_seq", True):
        is_xfmr = ppc.get("is_transformer")
        if is_xfmr is None:
            is_xfmr = np.zeros(ppc["branch"].shape[0], dtype=bool)
        p0["branch"][np.asarray(is_xfmr, dtype=bool), BR_STATUS] = 0
    return p0


def _gfm_shunt(ppc: dict[str, Any], gfm_buses: np.ndarray) -> np.ndarray:
    """Shunt admittance of grid-forming converters: 1 / Zvirtual, on system base.

    A GFM converter is a VOLTAGE source behind its virtual impedance, not a
    current injection. Electrically it therefore enters the fault network the
    same way a synchronous machine does -- and that is precisely why it restores
    a voltage reference to a network that has no synchronous plant left.

    The virtual impedance is on the converter's own MVA base, so it is referred
    to the system base the same way a machine reactance is.
    """
    n = ppc["bus"].shape[0]
    ysh = np.zeros(n, dtype=complex)
    zv_cfg = scenarios()["dynamics"]["gfm"]["virtual_impedance_pu"]
    gens = gen_at_bus(ppc)

    from .ibr import inverter_rating_pu

    for i in gfm_buses:
        rows = gens.get(int(i), [])
        if not rows:
            continue
        s_rating_pu = inverter_rating_pu(ppc, rows)      # already on system base
        if s_rating_pu <= 0:
            continue
        # converter base -> system base: Z_sys = Z_pu / S_rating_pu
        z_sys = complex(zv_cfg["r"], zv_cfg["x"]) / s_rating_pu
        if abs(z_sys) > 1e-12:
            ysh[i] += 1.0 / z_sys
    return ysh


def build_sequence_networks(
    ppc: dict[str, Any],
    v_prefault: np.ndarray,
    ibr_buses: np.ndarray | None = None,
    gfm_buses: np.ndarray | None = None,
) -> SequenceNetworks:
    """Assemble Y1/Y2/Y0 and their inverses for a fault study.

    Three source types, treated according to what they physically are:

      synchronous machines  voltage source behind X"d          -> shunt + Norton source
      GRID-FORMING (GFM)    voltage source behind Zvirtual     -> shunt + Norton source
      GRID-FOLLOWING (GFL)  current source, no internal EMF    -> no shunt; injected later

    GFL buses are deliberately given no shunt: an inverter that synchronises via
    a PLL has no internal voltage behind an impedance, and pretending it does
    would silently restore the very assumption this project tests.

    `gfm_buses` must be a subset of `ibr_buses`; anything in `ibr_buses` that is
    not grid-forming is treated as grid-following.
    """
    md = scenarios()["machine_data"]
    n = ppc["bus"].shape[0]
    ibr_buses = np.array([], dtype=int) if ibr_buses is None else np.asarray(ibr_buses, dtype=int)
    gfm_buses = np.array([], dtype=int) if gfm_buses is None else np.asarray(gfm_buses, dtype=int)

    ibr_set = set(ibr_buses.tolist())
    gfm_set = set(gfm_buses.tolist()) & ibr_set
    gfm_buses = np.array(sorted(gfm_set), dtype=int)
    gfl_buses = np.array(sorted(ibr_set - gfm_set), dtype=int)

    all_gen_buses = np.array(sorted(gen_at_bus(ppc).keys()), dtype=int)
    sync_buses = np.array([b for b in all_gen_buses if b not in ibr_set], dtype=int)

    y_load = _load_shunt(ppc, v_prefault)
    y_gfm = _gfm_shunt(ppc, gfm_buses)
    diag = np.arange(n)

    # --- positive sequence ---------------------------------------------------
    y1 = yb.build(ppc).astype(complex)
    ysh_1 = _machine_shunt(ppc, md["xdpp_pu"], md["r_over_x"], sync_buses)
    y1[diag, diag] += ysh_1 + y_gfm + y_load

    # --- negative sequence ---------------------------------------------------
    # A GFM converter presents its virtual impedance to negative sequence too,
    # but its control actively suppresses negative-sequence current, so the
    # effective impedance is higher. Modelled as 2x, a common approximation.
    y2 = yb.build(ppc).astype(complex)
    ysh_2 = _machine_shunt(ppc, md["x2_pu"], md["r_over_x"], sync_buses)
    y2[diag, diag] += ysh_2 + y_gfm / 2.0 + y_load

    # --- zero sequence -------------------------------------------------------
    # Converters are almost always interfaced through a delta-wye transformer,
    # so they contribute NO zero-sequence source. GFM gets no zero-seq shunt.
    y0 = yb.build(_zero_sequence_ppc(ppc)).astype(complex)
    x0_eff = md["x0_pu"] + 3.0 * md["xn_pu"]        # neutral impedance appears as 3*Zn
    ysh_0 = _machine_shunt(ppc, x0_eff, md["r_over_x"], sync_buses)
    y0[diag, diag] += ysh_0
    # a zero-sequence network with no ground path anywhere is singular
    y0[diag, diag] += 1e-9

    # Norton sources for every VOLTAGE source: I = E * Y, with E from prefault V
    e_int = np.zeros(n, dtype=complex)
    e_int[sync_buses] = v_prefault[sync_buses] * ysh_1[sync_buses]
    if gfm_buses.size:
        e_int[gfm_buses] = v_prefault[gfm_buses] * y_gfm[gfm_buses]

    return SequenceNetworks(
        y1=y1, y2=y2, y0=y0,
        z1=np.linalg.inv(y1), z2=np.linalg.inv(y2), z0=np.linalg.inv(y0),
        e_internal=e_int, ibr_buses=ibr_buses, sync_buses=sync_buses,
        v_prefault=v_prefault, gfm_buses=gfm_buses, gfl_buses=gfl_buses,
    )


# =============================================================================
# fault boundary conditions
# =============================================================================
def _fault_currents_from_oc(
    fault: FaultType,
    v1_oc: complex, v2_oc: complex, v0_oc: complex,
    z1: complex, z2: complex, z0: complex,
    zf: complex = 0.0, zg: complex = 0.0,
) -> tuple[complex, complex, complex]:
    """Sequence fault currents (I1, I2, I0) at the faulted bus.

    Generalised to nonzero negative/zero-sequence open-circuit voltage, which is
    what IBR control injection produces. Setting v2_oc = v0_oc = 0 recovers every
    textbook formula exactly -- `test_faults.py` asserts that.
    """
    za, zb, zc = z1 + zf, z2 + zf, z0 + zf + 3.0 * zg

    if fault == "3LG":
        # three phases tied together: V1' = V2' = 0, no ground path -> I0 = 0
        return v1_oc / za, v2_oc / zb, 0.0 + 0j

    if fault == "SLG":
        # phase a to ground: I1 = I2 = I0, V0'+V1'+V2' = 0
        i = (v1_oc + v2_oc + v0_oc) / (z1 + z2 + z0 + 3.0 * zf + 3.0 * zg)
        return i, i, i

    if fault == "LL":
        # phases b-c, no ground: I0 = 0, I2 = -I1, V1' = V2'
        i1 = (v1_oc - v2_oc) / (z1 + z2 + zf)
        return i1, -i1, 0.0 + 0j

    if fault == "LLG":
        # phases b-c to ground: I1+I2+I0 = 0, V1' = V2' = V0'
        mat = np.array([
            [1.0, 1.0, 1.0],
            [-za, zb, 0.0],
            [-za, 0.0, zc],
        ], dtype=complex)
        rhs = np.array([0.0, v2_oc - v1_oc, v0_oc - v1_oc], dtype=complex)
        i1, i2, i0 = np.linalg.solve(mat, rhs)
        return complex(i1), complex(i2), complex(i0)

    raise ValueError(f"unknown fault type {fault!r}")


def seq_to_phase(i1: complex, i2: complex, i0: complex) -> np.ndarray:
    """[I0, I1, I2] -> [Ia, Ib, Ic]."""
    return A_MATRIX @ np.array([i0, i1, i2], dtype=complex)


# =============================================================================
# results container
# =============================================================================
@dataclass
class FaultResult:
    fault_type: str
    bus: int
    method: str                      # 'classical' | 'ibr_aware'
    i1: complex = 0j
    i2: complex = 0j
    i0: complex = 0j
    iterations: int = 0
    converged: bool = True
    well_posed: bool = True
    note: str = ""
    v1: np.ndarray = field(default_factory=lambda: np.array([]))
    v2: np.ndarray = field(default_factory=lambda: np.array([]))
    v0: np.ndarray = field(default_factory=lambda: np.array([]))
    ibr_current: dict[int, complex] = field(default_factory=dict)

    @property
    def i_phase(self) -> np.ndarray:
        return seq_to_phase(self.i1, self.i2, self.i0)

    @property
    def i_fault_mag(self) -> float:
        """Magnitude of the largest faulted-phase current, per-unit."""
        return float(np.max(np.abs(self.i_phase)))

    @property
    def i_fault_angle_deg(self) -> float:
        """Angle of the largest faulted-phase current -- what directional relays see.

        Wrapped to (-180, 180] so that differences between two solutions are the
        true angular shift rather than a 360-degree artefact.
        """
        ph = self.i_phase
        from .metrics import wrap_deg
        return float(wrap_deg(np.degrees(np.angle(ph[int(np.argmax(np.abs(ph)))]))))

    @property
    def negative_sequence_ratio(self) -> float:
        """|I2|/|I1| -- the quantity negative-sequence relay elements depend on."""
        return float(abs(self.i2) / abs(self.i1)) if abs(self.i1) > 1e-12 else float("nan")


# =============================================================================
# classical solver
# =============================================================================
def solve_fault_classical(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    bus: int,
    fault: FaultType = "3LG",
    zf: complex = 0.0,
    zg: complex = 0.0,
) -> FaultResult:
    """Textbook symmetrical-component fault solve. Closed form, no iteration.

    Prefault voltage is taken from the load flow (rather than assuming a flat
    1.0 pu), which is the standard "superposition" formulation.
    """
    k = bus_index(ppc)[int(bus)]
    v1_oc = complex(seq.v_prefault[k])
    z1, z2, z0 = seq.z1[k, k], seq.z2[k, k], seq.z0[k, k]

    i1, i2, i0 = _fault_currents_from_oc(fault, v1_oc, 0j, 0j, z1, z2, z0, zf, zg)

    e_k = np.zeros(seq.n, dtype=complex)
    e_k[k] = 1.0
    return FaultResult(
        fault_type=fault, bus=int(bus), method="classical",
        i1=i1, i2=i2, i0=i0,
        v1=seq.v_prefault - seq.z1 @ (e_k * i1),
        v2=-(seq.z2 @ (e_k * i2)),
        v0=-(seq.z0 @ (e_k * i0)),
    )


# =============================================================================
# IBR-aware solver
# =============================================================================
def ibr_fault_injection(
    v1: complex, v2: complex,
    p_pu: float,
    i_limit: float,
    k2: float,
    reactive_priority: bool = True,
    v_prefault: complex | None = None,
) -> tuple[complex, complex]:
    """Current an IBR injects during a fault, per IEEE Std 2800-2022 in spirit.

    Positive sequence: reactive current support proportional to the voltage dip,
    Iq = K1 * (1 - |V1|), with active current taking whatever headroom is left
    under the converter limit (reactive priority) or the reverse.

    Negative sequence: I2 = -K2 * V2. This is the clause that matters most for
    protection -- without it (K2 = 0, i.e. pre-2800 inverters) the negative
    sequence current a distance/directional relay needs simply is not there.

    The hard clamp on |I| is what makes the network nonlinear, and is the entire
    reason the classical closed-form solution stops being valid.
    """
    vm = abs(v1)
    k1 = 2.0                                   # reactive support gain, pu/pu dip
    iq = min(k1 * max(0.0, 1.0 - vm), i_limit)  # capacitive support during a dip

    if reactive_priority:
        headroom = max(0.0, i_limit ** 2 - iq ** 2) ** 0.5
        ip = min(p_pu / vm if vm > 0.05 else 0.0, headroom)
    else:
        ip = min(p_pu / vm if vm > 0.05 else 0.0, i_limit)
        headroom = max(0.0, i_limit ** 2 - ip ** 2) ** 0.5
        iq = min(iq, headroom)

    # --- angle reference: PLL tracking, with coast-through on collapse -------
    # A grid-following converter injects current at an angle set by its PLL,
    # normally locked to the terminal voltage. During a close-in bolted fault
    # the terminal voltage can collapse to ~1e-5 pu, at which point V/|V| is
    # numerically meaningless and jitters with every iteration -- so no fixed
    # point exists and the solver cannot converge.
    #
    # Real inverters do not behave that way: below roughly 0.1 pu the PLL loses
    # lock and COASTS on its last angle. Modelling that coast is both physically
    # correct and what makes the problem well-conditioned. The two references
    # are blended over 0.02-0.10 pu so the injection stays continuous.
    v_hold = v_prefault if v_prefault is not None else (1.0 + 0j)
    ref_hold = v_hold / abs(v_hold) if abs(v_hold) > 1e-12 else 1.0 + 0j
    v_lo, v_hi = 0.02, 0.10

    if vm >= v_hi:
        ref = v1 / vm
    elif vm <= v_lo:
        ref = ref_hold                       # PLL has lost lock: coast
    else:
        w = (vm - v_lo) / (v_hi - v_lo)      # smooth handover
        blend = w * (v1 / vm) + (1.0 - w) * ref_hold
        ref = blend / abs(blend) if abs(blend) > 1e-12 else ref_hold

    i1 = ref * complex(ip, iq)

    # negative-sequence injection opposes V2 (IEEE 2800 unbalanced fault response)
    i2 = -k2 * v2
    if abs(i2) > i_limit:
        i2 = i2 / abs(i2) * i_limit
    return i1, i2


def solve_fault_ibr_aware(
    ppc: dict[str, Any],
    seq: SequenceNetworks,
    bus: int,
    fault: FaultType = "3LG",
    zf: complex = 0.0,
    zg: complex = 0.0,
    i_limit: float | None = None,
    k2: float | None = None,
    tol: float = 1e-8,
    max_iter: int = 200,
    relax: float = 0.7,
) -> FaultResult:
    """Iterative fault solve with current-limited IBRs.

    Fixed-point iteration under-relaxed by `relax`:
      1. From present IBR terminal voltages, evaluate the injected I1, I2.
      2. Superpose those injections onto the network to get open-circuit
         sequence voltages at the fault bus.
      3. Apply the fault boundary conditions -> new fault currents.
      4. Recompute bus voltages; repeat until the injections stop moving.

    Non-convergence is reported rather than hidden: the region where this loop
    fails to settle maps the boundary of validity of phasor fault analysis, and
    that boundary is a result worth plotting.
    """
    cfg = scenarios()["ibr"]
    i_limit = cfg["current_limit_pu"] if i_limit is None else i_limit
    k2 = cfg["ieee2800"]["k2"] if k2 is None else k2
    reactive_priority = cfg["ieee2800"].get("reactive_priority_during_fault", True)

    k = bus_index(ppc)[int(bus)]
    n = seq.n
    z1kk, z2kk, z0kk = seq.z1[k, k], seq.z2[k, k], seq.z0[k, k]

    # Well-posedness requires at least one VOLTAGE source. Synchronous machines
    # provide one; so do grid-forming converters. With neither, every source is
    # a current injection whose magnitude depends on a voltage that nothing
    # establishes, and the phasor fault problem is genuinely ill-posed rather
    # than merely hard to solve.
    #
    # This is the mechanism behind the WP7 mitigation result: converting part of
    # an all-grid-following fleet to grid-forming restores a voltage reference
    # and makes the problem solvable again.
    if seq.voltage_source_buses.size == 0:
        return FaultResult(
            fault_type=fault, bus=int(bus), method="ibr_aware",
            iterations=0, converged=False, well_posed=False,
            v1=seq.v_prefault.copy(),
            v2=np.zeros(n, dtype=complex), v0=np.zeros(n, dtype=complex),
            note="no voltage source remains (no synchronous plant, no grid-forming "
                 "converters) -- phasor fault model is ill-posed; add GFM or run EMT",
        )

    # Only GRID-FOLLOWING units iterate as nonlinear current injections.
    # Grid-forming units are already in the network as voltage sources behind
    # their virtual impedance, so they contribute through e_internal instead.
    gfl_buses = seq.gfl_buses if seq.gfl_buses.size or seq.gfm_buses.size else seq.ibr_buses

    gens = gen_at_bus(ppc)
    p_sched = {
        int(i): float(ppc["gen"][gens[int(i)], 1].sum()) / ppc["baseMVA"]
        for i in gfl_buses
    }

    # voltage-source contribution (synchronous + GFM) is linear and fixed
    v1_sync = seq.z1 @ seq.e_internal

    e_k = np.zeros(n, dtype=complex)
    e_k[k] = 1.0

    def evaluate(inj1: np.ndarray, inj2: np.ndarray):
        """Solve the network for a given set of GFL injections."""
        v1_oc_vec = v1_sync + seq.z1 @ inj1
        v2_oc_vec = seq.z2 @ inj2
        i1, i2, i0 = _fault_currents_from_oc(
            fault, complex(v1_oc_vec[k]), complex(v2_oc_vec[k]), 0j,
            z1kk, z2kk, z0kk, zf, zg,
        )
        v1 = v1_oc_vec - seq.z1 @ (e_k * i1)
        v2 = v2_oc_vec - seq.z2 @ (e_k * i2)
        v0 = -(seq.z0 @ (e_k * i0))
        return v1, v2, v0, i1, i2, i0

    def desired(v1: np.ndarray, v2: np.ndarray):
        """What each GFL converter WANTS to inject at these terminal voltages."""
        d1, d2 = np.zeros(n, dtype=complex), np.zeros(n, dtype=complex)
        for i in gfl_buses:
            a, b = ibr_fault_injection(
                v1[i], v2[i], p_sched.get(int(i), 0.0), i_limit, k2,
                reactive_priority, v_prefault=complex(seq.v_prefault[i]),
            )
            d1[i], d2[i] = a, b
        return d1, d2

    # ------------------------------------------------------------------
    # Unknowns are the GFL injections, packed as a real vector
    # x = [Re(I1), Im(I1), Re(I2), Im(I2)] over the GFL buses only.
    # There are at most a handful of GFL buses, so this vector is small and a
    # numerical Jacobian is cheap.
    # ------------------------------------------------------------------
    m = gfl_buses.size

    def pack(inj1: np.ndarray, inj2: np.ndarray) -> np.ndarray:
        if m == 0:
            return np.zeros(0)
        a, b = inj1[gfl_buses], inj2[gfl_buses]
        return np.concatenate([a.real, a.imag, b.real, b.imag])

    def unpack(x: np.ndarray):
        inj1, inj2 = np.zeros(n, dtype=complex), np.zeros(n, dtype=complex)
        if m:
            inj1[gfl_buses] = x[0:m] + 1j * x[m:2 * m]
            inj2[gfl_buses] = x[2 * m:3 * m] + 1j * x[3 * m:4 * m]
        return inj1, inj2

    def residual(x: np.ndarray) -> np.ndarray:
        """F(x) = desired_injection(V(x)) - x. A root of F is the solution."""
        inj1, inj2 = unpack(x)
        v1, v2, _, _, _, _ = evaluate(inj1, inj2)
        d1, d2 = desired(v1, v2)
        return pack(d1, d2) - x

    # With no grid-following units the problem is LINEAR: the grid-forming and
    # synchronous sources are already in the admittance matrix, so one solve
    # is exact and no iteration is needed at all.
    if m == 0:
        v1, v2, v0, i1, i2, i0 = evaluate(
            np.zeros(n, dtype=complex), np.zeros(n, dtype=complex))
        return FaultResult(
            fault_type=fault, bus=int(bus), method="ibr_aware",
            i1=i1, i2=i2, i0=i0, iterations=1, converged=True,
            v1=v1, v2=v2, v0=v0, ibr_current={},
            note="all voltage-source (GFM/synchronous): linear, solved directly",
        )

    x = np.zeros(4 * m)
    converged = False
    it = 0

    # --- Newton with a numerical Jacobian ------------------------------------
    # The fixed-point iteration this replaces stalled badly: the injection is a
    # hard-clamped function of the very voltage it produces, so the map is only
    # weakly contractive and roughly half of all fault cases hit the iteration
    # cap. Newton converges on the same problem in a handful of steps.
    f = residual(x)
    for it in range(1, max_iter + 1):
        if np.max(np.abs(f)) < tol:
            converged = True
            break

        eps = 1e-7
        jac = np.empty((4 * m, 4 * m))
        for j in range(4 * m):
            xp = x.copy()
            xp[j] += eps
            jac[:, j] = (residual(xp) - f) / eps

        try:
            dx = np.linalg.solve(jac, -f)
        except np.linalg.LinAlgError:
            dx = -f * relax                      # fall back to a damped step

        # backtracking line search keeps the clamped nonlinearity from
        # throwing Newton past the solution
        step, f_new = 1.0, None
        for _ in range(20):
            cand = residual(x + step * dx)
            if np.max(np.abs(cand)) < np.max(np.abs(f)):
                f_new = cand
                break
            step *= 0.5
        if f_new is None:                        # no improvement: damped update
            step = relax
            f_new = residual(x + step * dx)

        x = x + step * dx
        f = f_new
    else:
        converged = bool(np.max(np.abs(f)) < tol)

    inj1, inj2 = unpack(x)
    v1, v2, v0, i1, i2, i0 = evaluate(inj1, inj2)

    return FaultResult(
        fault_type=fault, bus=int(bus), method="ibr_aware",
        i1=i1, i2=i2, i0=i0, iterations=it, converged=converged,
        v1=v1, v2=v2, v0=v0,
        ibr_current={int(i): complex(inj1[i]) for i in gfl_buses},
    )


# =============================================================================
# short-circuit MVA (feeds the system-strength module)
# =============================================================================
def short_circuit_mva(
    ppc: dict[str, Any], seq: SequenceNetworks, bus: int,
) -> float:
    """Three-phase short-circuit level at a bus, in MVA.

    Ssc = |V_pre|^2 / |Z1_thevenin| on a per-unit basis, times the system base.
    This is the numerator of every system-strength metric in `strength.py`.
    """
    k = bus_index(ppc)[int(bus)]
    z_th = seq.z1[k, k]
    if abs(z_th) < 1e-12:
        return float("inf")
    return float(abs(seq.v_prefault[k]) ** 2 / abs(z_th) * ppc["baseMVA"])
