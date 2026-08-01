"""Regression gates for the IEEE IBR benchmark.

These are not unit tests for their own sake -- each one guards a claim the report
makes. If a gate fails, the corresponding result in the report is no longer
supported and must not be published.

    cd IEEE/python
    python -m pytest tests/ -v
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import (
    config, estimation, faults as F, ibr as I, metrics, ppc as P, solvers as S,
    stability, strength, transient,
)

CASES = ["case9", "case14", "case30", "case39"]


# =============================================================================
# foundation
# =============================================================================
@pytest.mark.parametrize("case", CASES)
def test_ybus_matches_pandapower(case):
    """Our from-scratch Ybus must equal pandapower's to machine precision.

    Everything downstream is built on Ybus. If this fails, nothing else matters.
    """
    from gridbench import ybus as Y
    err = Y.validate_against_reference(P.load_ppc(case), atol=1e-9)
    assert err < 1e-12, f"Ybus mismatch {err:.3e}"


@pytest.mark.parametrize("case", CASES)
def test_bus_numbering_is_one_based(case):
    c = P.load_ppc(case)
    buses = c["bus"][:, P.BUS_I].astype(int)
    assert buses.min() == 1
    assert set(buses) == set(range(1, len(buses) + 1))


# =============================================================================
# WP1 -- solvers
# =============================================================================
@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("alg", ["nr", "fdxb", "fdbx"])
def test_solver_matches_pandapower(case, alg):
    """Hand-written NR/FDLF must reproduce pandapower to 1e-8 pu."""
    import pandapower as pp
    import pandapower.networks as pn

    net = getattr(pn, case)()
    pp.runpp(net, numba=False)

    r = S.solve(P.load_ppc(case), alg)
    assert r.converged, f"{alg} diverged on {case}"

    tol_vm = config.tol("cross_tool", "handwritten_vs_matpower", "vm_pu")
    assert metrics.max_abs_error(r.vm, net.res_bus.vm_pu.values) <= tol_vm


@pytest.mark.parametrize("case", CASES)
def test_newton_converges_quadratically(case):
    """NR should reach 1e-10 in a handful of iterations regardless of size."""
    r = S.newton_raphson(P.load_ppc(case))
    assert r.converged
    assert r.iterations <= 6, f"{case} took {r.iterations} iterations"


def test_gauss_seidel_is_acceleration_sensitive():
    """Documents the WP1 finding: GS on case39 depends on a parameter NR lacks."""
    c = P.load_ppc("case39")
    assert S.gauss_seidel(c, accel=1.4, max_iter=8000).converged
    assert not S.gauss_seidel(c, accel=1.8, max_iter=8000).converged


# =============================================================================
# WP2 -- the experimental control
# =============================================================================
@pytest.mark.parametrize("case", ["case14", "case30", "case39"])
def test_inverter_sizing_matches_machine_capability(case):
    """THE control for WP2.

    At nominal voltage and rated power the inverter's reactive capability must
    equal the displaced machine's case-file QMAX. If it does not, the classical
    vs IBR-aware comparison measures the sizing assumption, not the physics.
    """
    c = P.load_ppc(case)
    for g in range(c["gen"].shape[0]):
        if c["is_slack_gen"][g] or c["gen"][g, P.GEN_STATUS] <= 0:
            continue
        s_rating = I.inverter_rating_pu(c, [g])
        pmax_pu = float(c["gen"][g, P.PMAX]) / c["baseMVA"]
        cap = I.q_capability_pu(1.0, pmax_pu, 1.0, s_rating) * c["baseMVA"]
        assert abs(cap - float(c["gen"][g, P.QMAX])) < 1e-6


@pytest.mark.parametrize("case", ["case14", "case30", "case39"])
def test_zero_penetration_reduces_to_classical(case):
    """At 0% IBR the two models must be the SAME computation, bit for bit."""
    c = P.load_ppc(case)
    classical = S.newton_raphson(c, tol=1e-11, enforce_q_limits=True)
    r = I.ibr_powerflow(c, np.array([], dtype=int))
    assert r.converged and classical.converged
    assert metrics.max_abs_error(r.vm, classical.vm) < 1e-12


def test_penetration_selection_respects_target():
    """Regression: the greedy selector once absorbed every zero-MW machine."""
    c = P.load_ppc("case14")
    gens = I.select_ibr_gens(c, 20.0)
    assert 0 < len(gens) < c["gen"].shape[0]
    assert I.actual_penetration_pct(c, gens) < 45.0


# =============================================================================
# WP3 -- faults
# =============================================================================
def test_classical_fault_matches_closed_form():
    """Zbus solver vs the textbook closed forms, on IEEE 14-bus."""
    c = P.load_ppc("case14")
    seq = F.build_sequence_networks(c, S.newton_raphson(c).v)
    k = P.bus_index(c)[4]
    z1, z2, z0 = seq.z1[k, k], seq.z2[k, k], seq.z0[k, k]
    vf = seq.v_prefault[k]

    checks = {
        "3LG": abs(vf / z1),
        "SLG": abs(3 * vf / (z1 + z2 + z0)),
        "LL": abs(np.sqrt(3) * vf / (z1 + z2)),
    }
    for ft, expected in checks.items():
        got = F.solve_fault_classical(c, seq, 4, ft).i_fault_mag
        assert abs(got - expected) / expected < 1e-9, f"{ft}: {got} vs {expected}"


def test_fault_type_ordering_is_physical():
    """3LG must be the most severe; LL should be ~0.866 of it."""
    c = P.load_ppc("case14")
    seq = F.build_sequence_networks(c, S.newton_raphson(c).v)
    mags = {ft: F.solve_fault_classical(c, seq, 4, ft).i_fault_mag
            for ft in ["3LG", "SLG", "LL", "LLG"]}
    assert mags["3LG"] > mags["LL"] > mags["SLG"]
    assert abs(mags["LL"] / mags["3LG"] - np.sqrt(3) / 2) < 0.02


def test_bolted_three_phase_fault_collapses_voltage():
    c = P.load_ppc("case14")
    seq = F.build_sequence_networks(c, S.newton_raphson(c).v)
    r = F.solve_fault_classical(c, seq, 4, "3LG")
    assert abs(r.v1[P.bus_index(c)[4]]) < 1e-9


def test_ibr_reduces_fault_current():
    """The core WP3 claim: current-limited inverters cut fault current."""
    c = P.load_ppc("case39")
    base = S.newton_raphson(c)
    gens = I.select_ibr_gens(c, 60.0)
    seq = F.build_sequence_networks(c, base.v, ibr_buses=I.ibr_buses_from_gens(c, gens))

    a = F.solve_fault_classical(c, seq, 16, "3LG")
    b = F.solve_fault_ibr_aware(c, seq, 16, "3LG")
    assert b.converged
    assert b.i_fault_mag < a.i_fault_mag * 0.9


def test_no_voltage_source_is_reported_ill_posed():
    """100% grid-FOLLOWING IBR has no voltage reference: flagged, not solved."""
    c = P.load_ppc("case39")
    base = S.newton_raphson(c)
    gens = I.select_ibr_gens(c, 100.0)
    seq = F.build_sequence_networks(c, base.v, ibr_buses=I.ibr_buses_from_gens(c, gens))
    assert seq.voltage_source_buses.size == 0
    r = F.solve_fault_ibr_aware(c, seq, 16, "3LG")
    assert not r.well_posed
    assert not r.converged
    assert "ill-posed" in r.note


# =============================================================================
# WP5 -- mitigation
# =============================================================================
def test_grid_forming_restores_well_posedness():
    """THE mitigation result: GFM converters are voltage sources, so they make
    a 100%-inverter network solvable again."""
    c = P.load_ppc("case39")
    base = S.newton_raphson(c)
    gens = I.select_ibr_gens(c, 100.0)
    ibr_buses = I.ibr_buses_from_gens(c, gens)

    # all grid-following -> ill-posed
    seq_gfl = F.build_sequence_networks(c, base.v, ibr_buses=ibr_buses)
    assert not F.solve_fault_ibr_aware(c, seq_gfl, 16, "3LG").well_posed

    # all grid-forming -> well-posed, and LINEAR (no current sources to iterate)
    _, gfm = I.split_gfl_gfm(c, gens, 100.0)
    seq_gfm = F.build_sequence_networks(
        c, base.v, ibr_buses=ibr_buses, gfm_buses=I.ibr_buses_from_gens(c, gfm))
    r = F.solve_fault_ibr_aware(c, seq_gfm, 16, "3LG")
    assert r.well_posed and r.converged
    assert r.i_fault_mag > 0


def test_grid_forming_share_monotonically_raises_strength():
    """More grid-forming capacity must mean more system strength, not less."""
    c = P.load_ppc("case39")
    base = S.newton_raphson(c)
    gens = I.select_ibr_gens(c, 100.0)
    ibr_buses = I.ibr_buses_from_gens(c, gens)

    prev = -np.inf
    for share in [0, 25, 50, 75, 100]:
        _, gfm = I.split_gfl_gfm(c, gens, share)
        seq = F.build_sequence_networks(
            c, base.v, ibr_buses=ibr_buses,
            gfm_buses=I.ibr_buses_from_gens(c, gfm) if gfm.size else None)
        wscr = strength.profile(c, seq, gens)["wscr"]
        assert wscr >= prev - 1e-9, f"WSCR fell at {share}% GFM"
        prev = wscr


def test_gfl_gfm_split_partitions_the_fleet():
    c = P.load_ppc("case39")
    gens = I.select_ibr_gens(c, 80.0)
    for share in [0, 25, 50, 75, 100]:
        gfl, gfm = I.split_gfl_gfm(c, gens, share)
        assert set(gfl.tolist()) | set(gfm.tolist()) == set(gens.tolist())
        assert not (set(gfl.tolist()) & set(gfm.tolist()))


def test_synchronous_condenser_raises_short_circuit_level():
    """A condenser must add fault current at the bus it is installed at."""
    c0 = P.load_ppc("case39")
    base0 = S.newton_raphson(c0)
    gens = I.select_ibr_gens(c0, 80.0)
    seq0 = F.build_sequence_networks(
        c0, base0.v, ibr_buses=I.ibr_buses_from_gens(c0, gens))

    gen_buses = {int(b) for b in c0["gen"][:, P.GEN_BUS]}
    target = min((int(b) for b in c0["bus"][:, P.BUS_I].astype(int)
                  if b not in gen_buses),
                 key=lambda b: F.short_circuit_mva(c0, seq0, b))
    before = F.short_circuit_mva(c0, seq0, target)

    c = I.add_synchronous_condensers(c0, [target], 300.0)
    base = S.newton_raphson(c, enforce_q_limits=True)
    seq = F.build_sequence_networks(
        c, base.v, ibr_buses=I.ibr_buses_from_gens(c, gens))
    after = F.short_circuit_mva(c, seq, target)

    assert after > before * 1.1, f"condenser had no effect: {before} -> {after}"


def test_condenser_is_not_selected_as_an_ibr():
    """Regression: condensers were once picked up by the IBR selector."""
    c0 = P.load_ppc("case39")
    n_before = c0["gen"].shape[0]
    c = I.add_synchronous_condensers(c0, [12], 200.0)
    gens = I.select_ibr_gens(c, 80.0)
    assert not (set(gens.tolist()) & set(range(n_before, c["gen"].shape[0])))


# =============================================================================
# metrics
# =============================================================================
def test_angle_wrapping():
    """A -161.4 -> 177.0 shift is -21.6 deg, not +338.5."""
    assert abs(metrics.angle_diff_deg(177.04, -161.43) - (-21.53)) < 0.01
    assert metrics.wrap_deg(360.0) == pytest.approx(0.0)
    assert metrics.wrap_deg(-180.0) == pytest.approx(180.0)


def test_pct_error_guards_zero_reference():
    assert np.isnan(metrics.pct_error(1.0, 0.0))


# =============================================================================
# WP5 / WP6
# =============================================================================
def test_inertia_falls_and_rocof_rises_with_penetration():
    c = P.load_ppc("case39")
    dp = ST_dp = 1000.0
    prev_h, prev_rocof = np.inf, 0.0
    for pen in [0, 20, 40, 60, 80]:
        gens = I.select_ibr_gens(c, pen)
        from gridbench import inertia
        fr = inertia.frequency_response(c, dp, gens)
        assert fr.h_sys_s <= prev_h + 1e-9
        assert abs(fr.rocof_hz_s) >= abs(prev_rocof) - 1e-9
        prev_h, prev_rocof = fr.h_sys_s, fr.rocof_hz_s


def test_n1_screen_separates_islanding_from_divergence():
    c = P.load_ppc("case39")
    results = stability.n1_screen(c)
    assert any(r.islanded for r in results), "case39 has radial generator branches"
    for r in results:
        assert not (r.islanded and r.converged)


def test_scr_falls_with_penetration():
    c = P.load_ppc("case39")
    base = S.newton_raphson(c)
    prev = np.inf
    for pen in [20, 40, 60, 80]:
        gens = I.select_ibr_gens(c, pen)
        seq = F.build_sequence_networks(
            c, base.v, ibr_buses=I.ibr_buses_from_gens(c, gens))
        wscr = strength.profile(c, seq, gens)["wscr"]
        assert wscr <= prev + 1e-9
        prev = wscr


# =============================================================================
# WP7 -- transient stability
# =============================================================================
def test_dynamic_network_is_at_equilibrium():
    """Pe must equal Pm at t=0, or the machines swing from a false start."""
    c = P.load_ppc("case9")
    pf = S.newton_raphson(c)
    net = transient.build_dynamic_network(c, pf.v)
    pe = transient._electrical_power(net.e_internal, net.delta0, net.y_reduced)
    assert np.allclose(pe, net.p_mech, atol=1e-6), "not at equilibrium"


def test_fault_reduces_power_transfer():
    """A bolted fault must reduce TOTAL power transfer and collapse the nearest unit.

    Not every machine individually: power redistributes during a fault, so a
    machine electrically distant from it can legitimately pick up output (on
    case9, a fault at bus 7 raises gen 1 from 0.720 to 0.759 pu). What must hold
    is that the system as a whole cannot deliver its pre-fault power, and that
    the unit nearest the fault is severely curtailed -- that mismatch is the
    accelerating power that drives the swing.
    """
    c = P.load_ppc("case9")
    pf = S.newton_raphson(c)
    net = transient.build_dynamic_network(c, pf.v)
    idx = P.bus_index(c)

    for fault_bus in (4, 6, 7, 8):
        y_f = transient.reduced_for_condition(
            c, pf.v, net, None, None, fault_bus=idx[fault_bus])
        pe = transient._electrical_power(net.e_internal, net.delta0, y_f)
        assert pe.sum() < net.p_mech.sum(), \
            f"fault at {fault_bus} did not reduce total transfer"
        worst = float(np.min(pe / net.p_mech))
        assert worst < 0.6, \
            f"fault at {fault_bus} left every unit above 60% output ({worst:.2f})"


def test_published_machine_data_is_used():
    """case9/case39 must use benchmark H, not the generic default."""
    c = P.load_ppc("case9")
    pf = S.newton_raphson(c)
    net = transient.build_dynamic_network(c, pf.v)
    # Anderson & Fouad, referred to 100 MVA: 58.5 / 12.29 / 3.85
    assert abs(net.h.sum() - 74.64) < 0.5, f"unexpected total inertia {net.h.sum()}"


def test_cct_is_bracketed_on_benchmark():
    """IEEE 39-bus CCT must be finite and in a physically sane range."""
    c = P.load_ppc("case39")
    pf = S.newton_raphson(c)
    br = c["branch"]
    k = next(j for j in range(br.shape[0])
             if {int(br[j, 0]), int(br[j, 1])} == {16, 17})
    net = transient.build_dynamic_network(c, pf.v)
    net.damping[:] = 0.0
    r = transient.critical_clearing_time(c, pf.v, net, 16, outaged_branch=k, t_max=1.0)
    assert r["bracketed"]
    assert 0.05 < r["cct_s"] < 0.5, f"CCT {r['cct_s']} outside plausible range"


def test_virtual_inertia_extends_cct():
    """THE controlled GFM result: more virtual inertia -> longer CCT."""
    c = P.load_ppc("case39")
    pf = S.newton_raphson(c)
    br = c["branch"]
    k = next(j for j in range(br.shape[0])
             if {int(br[j, 0]), int(br[j, 1])} == {16, 17})
    gens = I.select_ibr_gens(c, 60.0)
    _, gfm = I.split_gfl_gfm(c, gens, 100.0)
    rows = transient.virtual_inertia_sweep(
        c, pf.v, 16, gens, gfm, [0.5, 4.0, 16.0], outaged_branch=k, t_max=1.0)

    # unit count must be constant -- that is what makes this controlled
    assert len({r["n_dynamic_units"] for r in rows}) == 1
    assert rows[-1]["cct_s"] > rows[0]["cct_s"] * 1.5, \
        f"virtual inertia did not extend CCT: {[r['cct_s'] for r in rows]}"


# =============================================================================
# WP8 -- state estimation and FDIA
# =============================================================================
def test_estimator_recovers_true_state():
    c = P.load_ppc("case14")
    pf = S.newton_raphson(c)
    meas = estimation.build_measurements(c, pf.v)
    r = estimation.estimate(c, meas)
    assert r.converged
    assert np.max(np.abs(np.abs(r.v) - np.abs(pf.v))) < 5e-3
    assert not r.bad_data_detected, "false alarm on clean data"


def test_random_attack_is_detected():
    """The defence must work against a naive attacker."""
    c = P.load_ppc("case14")
    pf = S.newton_raphson(c)
    meas = estimation.build_measurements(c, pf.v)
    a = estimation.random_attack(meas, n_targets=3, magnitude=0.2)
    r = estimation.estimate(c, estimation.apply_attack(meas, a))
    assert r.bad_data_detected


def test_exact_fdia_is_undetectable():
    """The core security result: an exact AC FDIA leaves the residual unchanged.

    Not 'hard to detect' -- structurally invisible to a residual test, at any
    magnitude. This is what motivates out-of-band verification.
    """
    c = P.load_ppc("case14")
    pf = S.newton_raphson(c)
    meas = estimation.build_measurements(c, pf.v)
    clean = estimation.estimate(c, meas)
    i = P.bus_index(c)[5]

    for deg in (5.0, 10.0, 20.0):
        a = estimation.stealthy_attack(c, meas, pf.v, 5, np.deg2rad(deg))
        r = estimation.estimate(c, estimation.apply_attack(meas, a))
        assert not r.bad_data_detected, f"{deg} deg attack was detected"
        # residual essentially unchanged from clean
        assert abs(r.objective - clean.objective) < 1.0
        # but the operator's state IS corrupted, by the intended amount
        err = abs(np.degrees(np.angle(r.v[i]) - np.angle(pf.v[i])))
        assert err > deg * 0.8, f"attack did not move the state ({err:.2f} deg)"


def test_linearised_fdia_is_detectable_at_scale():
    """The DC-model attack is NOT stealthy against an AC estimator."""
    c = P.load_ppc("case39")
    pf = S.newton_raphson(c)
    meas = estimation.build_measurements(c, pf.v)
    a = estimation.stealthy_attack(c, meas, pf.v, 8, np.deg2rad(10.0), linear=True)
    r = estimation.estimate(c, estimation.apply_attack(meas, a))
    assert r.bad_data_detected, "linearised attack should leave a residual"
