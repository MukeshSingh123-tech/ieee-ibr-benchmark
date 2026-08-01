"""Low-inertia frequency analysis: system inertia, RoCoF, nadir, and the
minimum-inertia constraint.

Inertia never used to be a planning variable. It came free with synchronous
generation, in proportion to it. Displace those machines with inverters and it
disappears -- an inverter has no rotating mass coupled to system frequency, so
100 MW of PV contributes exactly zero to the inertial response that 100 MW of
steam plant provided.

The consequences are governed by the swing equation for the system centre of
inertia. For a sudden generation loss dP at t = 0:

        2 H_sys  d(df)/dt = -dP           =>   RoCoF_0 = -dP * f0 / (2 H_sys)

RoCoF is therefore INVERSELY proportional to system inertia: halve the inertia
and the initial rate of change of frequency doubles. Since RoCoF-based
loss-of-mains protection trips at a fixed threshold (typically 0.5-1.0 Hz/s),
there is a hard penetration limit beyond which a credible contingency trips
protection across the system -- which is a cascading-failure mechanism, not a
power-quality inconvenience.

Everything here is the standard low-order System Frequency Response model. It is
a reduced model, not an EMT simulation, and its assumptions are stated in
`docs/03_methodology.md`: single centre-of-inertia frequency, aggregated
governor, no network dynamics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .config import scenarios
from .ppc import GEN_BUS, GEN_STATUS, PG, PMAX


# =============================================================================
# system inertia
# =============================================================================
def machine_inertia_mws(ppc: dict[str, Any], gen_row: int, h_seconds: float) -> float:
    """Stored kinetic energy of one machine, in MW-s.

    E = H * S_rating, where H is the inertia constant in seconds (energy stored
    at rated speed divided by rating). Typical H: 2-4 s for hydro, 3-6 s for
    steam, 4-8 s for large thermal units.
    """
    rating = ppc["gen_mva_base"][gen_row]
    if not np.isfinite(rating) or rating <= 0:
        rating = abs(float(ppc["gen"][gen_row, PG])) or ppc["baseMVA"]
    return float(h_seconds * rating)


def system_inertia(
    ppc: dict[str, Any],
    ibr_gens: np.ndarray | None = None,
    h_seconds: float | None = None,
) -> dict[str, float]:
    """Aggregate system inertia with IBR-displaced machines contributing nothing.

    Returns both the stored energy (MW-s) and the normalised H_sys in seconds,
    referred to total online generation -- the form the swing equation needs.
    """
    cfg = scenarios().get("dynamics", {})
    h_seconds = float(cfg.get("h_default_s", 4.0)) if h_seconds is None else h_seconds

    gen = ppc["gen"]
    live = np.flatnonzero(gen[:, GEN_STATUS] > 0)
    ibr = set(np.asarray(ibr_gens, dtype=int).tolist()) if ibr_gens is not None else set()

    energy_mws = 0.0
    sync_mva = 0.0
    for g in live:
        if int(g) in ibr:
            continue                      # an inverter stores no rotational energy
        energy_mws += machine_inertia_mws(ppc, int(g), h_seconds)
        rating = ppc["gen_mva_base"][int(g)]
        sync_mva += float(rating) if np.isfinite(rating) and rating > 0 else 0.0

    total_gen_mw = float(gen[live, PG].sum())
    h_sys = energy_mws / total_gen_mw if total_gen_mw > 0 else 0.0
    return {
        "energy_mws": energy_mws,
        "h_sys_s": h_sys,
        "synchronous_mva": sync_mva,
        "total_gen_mw": total_gen_mw,
    }


# =============================================================================
# frequency response
# =============================================================================
@dataclass
class FrequencyResponse:
    rocof_hz_s: float
    nadir_hz: float
    nadir_time_s: float
    settling_hz: float
    h_sys_s: float
    delta_p_mw: float
    violates_rocof: bool
    violates_nadir: bool

    def summary(self) -> str:
        return (
            f"H={self.h_sys_s:.2f}s  RoCoF={self.rocof_hz_s:.3f} Hz/s  "
            f"nadir={self.nadir_hz:.3f} Hz @ {self.nadir_time_s:.2f}s"
        )


def initial_rocof(delta_p_mw: float, h_sys_s: float, total_gen_mw: float,
                  f0: float = 60.0) -> float:
    """Initial rate of change of frequency from the swing equation, Hz/s.

    RoCoF = -dP * f0 / (2 * H_sys * S_base). Negative for a generation loss.
    """
    denom = 2.0 * h_sys_s * total_gen_mw
    if denom <= 0:
        return float("-inf")
    return float(-delta_p_mw * f0 / denom)


def frequency_response(
    ppc: dict[str, Any],
    delta_p_mw: float,
    ibr_gens: np.ndarray | None = None,
    h_seconds: float | None = None,
) -> FrequencyResponse:
    """Low-order System Frequency Response (SFR) model of a generation-loss event.

    Uses the standard Anderson-Mirheydar closed form for a reheat-turbine system:
    a second-order response whose nadir depends on inertia, governor droop, and
    reheat time constant. The nadir -- not RoCoF -- is what determines whether
    under-frequency load shedding operates.
    """
    cfg = scenarios().get("dynamics", {})
    f0 = float(scenarios()["meta"].get("base_frequency_hz", 60.0))
    r = float(cfg.get("governor_droop_pu", 0.05))       # 5% droop
    tr = float(cfg.get("reheat_time_s", 8.0))
    km = float(cfg.get("mechanical_power_gain", 0.95))
    fh = float(cfg.get("high_pressure_fraction", 0.3))

    d_load = float(cfg.get("load_damping_pu", 1.0))

    inertia = system_inertia(ppc, ibr_gens, h_seconds)
    h = inertia["h_sys_s"]
    s_base = inertia["total_gen_mw"]

    rocof = initial_rocof(delta_p_mw, h, s_base, f0)
    pu_step = delta_p_mw / s_base if s_base > 0 else 0.0

    if h <= 0 or s_base <= 0:
        # no synchronous inertia at all: RoCoF is unbounded and the phasor SFR
        # model has nothing left to integrate. Reported, not silently smoothed.
        return FrequencyResponse(float("-inf"), float("-inf"), float("nan"),
                                 float("-inf"), h, delta_p_mw, True, True)

    # Integrate the SFR model directly rather than using the Anderson-Mirheydar
    # closed form. That closed form assumes a moderately damped second-order
    # response and returns nonsense as H approaches zero -- it reported the nadir
    # IMPROVING at 84% penetration, which is backwards. Numerical integration
    # stays valid across the whole penetration sweep, which is the entire region
    # of interest here.
    #
    #   2H d(df)/dt = -dP + dPm - D*df
    #   dPm         = -(Km/R) * [ Fh*df + (1-Fh)*x ]
    #   Tr dx/dt    = df - x                      (reheat lag)
    def derivatives(state: np.ndarray) -> np.ndarray:
        df, x = state
        dpm = -(km / r) * (fh * df + (1.0 - fh) * x)
        ddf = (-pu_step + dpm - d_load * df) / (2.0 * h)
        dx = (df - x) / tr
        return np.array([ddf, dx])

    dt = min(0.01, tr / 200.0, h / 10.0)
    n_steps = int(60.0 / dt)                 # 60 s is well past the nadir
    state = np.zeros(2)
    nadir_pu, t_nadir = 0.0, 0.0

    for k in range(n_steps):
        k1 = derivatives(state)
        k2 = derivatives(state + 0.5 * dt * k1)
        k3 = derivatives(state + 0.5 * dt * k2)
        k4 = derivatives(state + dt * k3)
        state = state + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if state[0] < nadir_pu:
            nadir_pu, t_nadir = float(state[0]), (k + 1) * dt
        if not np.isfinite(state).all():
            break

    settling_pu = float(state[0])
    nadir_hz = f0 + nadir_pu * f0
    settling_hz = f0 + settling_pu * f0

    limits = scenarios().get("frequency_limits", {})
    rocof_limit = float(limits.get("rocof_hz_s", 1.0))
    ufls_hz = float(limits.get("ufls_threshold_hz", 59.3))

    return FrequencyResponse(
        rocof_hz_s=rocof,
        nadir_hz=float(nadir_hz),
        nadir_time_s=float(t_nadir),
        settling_hz=float(settling_hz),
        h_sys_s=h,
        delta_p_mw=delta_p_mw,
        violates_rocof=abs(rocof) > rocof_limit,
        violates_nadir=nadir_hz < ufls_hz,
    )


def minimum_inertia_mws(
    delta_p_mw: float, rocof_limit_hz_s: float, f0: float = 60.0,
) -> float:
    """Stored energy required to keep RoCoF within limit for a given contingency.

        E_min = dP * f0 / (2 * RoCoF_limit)

    This is the constraint system operators now enforce in dispatch (EirGrid,
    AEMO, National Grid all have a version of it). Below E_min the operator must
    either commit synchronous plant, add synchronous condensers, procure
    fast-frequency response, or curtail the interconnector -- all of which cost
    money, which is why the number matters commercially as well as technically.
    """
    if rocof_limit_hz_s <= 0:
        return float("inf")
    return float(abs(delta_p_mw) * f0 / (2.0 * rocof_limit_hz_s))


def max_ibr_penetration_for_rocof(
    ppc: dict[str, Any],
    delta_p_mw: float,
    rocof_limit_hz_s: float | None = None,
    h_seconds: float | None = None,
) -> float:
    """Highest IBR penetration (%) that still respects the RoCoF limit.

    Found by bisection over the penetration sweep, using the same displacement
    ordering as the rest of the project so the answer is comparable with the
    power-flow and fault results at the same penetration levels.
    """
    from .ibr import actual_penetration_pct, select_ibr_gens

    limits = scenarios().get("frequency_limits", {})
    rocof_limit = float(limits.get("rocof_hz_s", 1.0)) if rocof_limit_hz_s is None \
        else rocof_limit_hz_s

    best = 0.0
    for pen in np.arange(0.0, 100.5, 2.5):
        gens = select_ibr_gens(ppc, float(pen))
        fr = frequency_response(ppc, delta_p_mw, gens, h_seconds)
        if abs(fr.rocof_hz_s) <= rocof_limit:
            best = actual_penetration_pct(ppc, gens)
        else:
            break
    return float(best)


def largest_single_contingency_mw(ppc: dict[str, Any], exclude_slack: bool = True) -> float:
    """Largest online unit -- the standard N-1 generation-loss contingency.

    The slack is excluded by default: in these cases it stands for the rest of
    the interconnection rather than a physical machine, so "losing" it is not a
    credible single contingency. On IEEE 14-bus including it would mean losing
    232 of 272 MW (85% of generation) and every frequency result would be
    dominated by that artefact rather than by inertia.
    """
    gen = ppc["gen"]
    live = gen[:, GEN_STATUS] > 0
    if exclude_slack:
        live = live & ~ppc["is_slack_gen"]
    return float(gen[live, PG].max()) if live.any() else 0.0
