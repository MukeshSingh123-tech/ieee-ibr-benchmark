"""WP3/WP4 -- fault analysis with current-limited IBRs, and the protection consequences.

Three questions, in order:

  WP3a  How far does classical symmetrical-component fault analysis deviate once
        the sources are current-limited inverters instead of synchronous machines?

  WP3b  Does IEEE Std 2800-2022 negative-sequence injection (I2 = K2*V2) change
        the answer? K2 = 0 represents a pre-2800 inverter, which injects no
        negative-sequence current at all.

  WP4   What do standard, correctly-set relays DO with those currents? The
        relays are commissioned from a classical study; the fault is then solved
        with the IBR-aware model, and the two verdicts are compared.

Outputs
    results/tables/wp3_fault_error.csv
    results/tables/wp3_ieee2800_k2_sweep.csv
    results/tables/wp4_protection_misoperation.csv
"""

from __future__ import annotations

import sys
import warnings
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import (  # noqa: E402
    config, faults as F, ibr as I, metrics, ppc as P, protection as PR, solvers as S,
)

CASES = ["case14", "case39"]
FAULT_TYPES = ["3LG", "SLG", "LL", "LLG"]


def _setup(case: str, pen: float):
    c = P.load_ppc(case)
    base = S.newton_raphson(c)
    gens = I.select_ibr_gens(c, pen)
    seq = F.build_sequence_networks(
        c, base.v, ibr_buses=I.ibr_buses_from_gens(c, gens),
    )
    return c, base, gens, seq


def fault_error_sweep() -> list[dict]:
    """Classical vs IBR-aware fault current, per fault type and penetration."""
    levels = config.scenarios()["penetration"]["levels_pct"]
    rows = []

    for case in CASES:
        locations = config.scenarios()["faults"]["locations"].get(case, [])
        for pen in levels:
            c, base, gens, seq = _setup(case, pen)
            if gens.size == 0:
                continue
            for bus in locations:
                for ft in FAULT_TYPES:
                    a = F.solve_fault_classical(c, seq, bus, ft)
                    b = F.solve_fault_ibr_aware(c, seq, bus, ft)
                    rows.append({
                        "case": case,
                        "pen_pct": I.actual_penetration_pct(c, gens),
                        "bus": bus,
                        "fault_type": ft,
                        "if_classical_pu": a.i_fault_mag,
                        "if_ibr_pu": b.i_fault_mag,
                        "error_pct": metrics.pct_error(b.i_fault_mag, a.i_fault_mag),
                        "angle_classical_deg": a.i_fault_angle_deg,
                        "angle_ibr_deg": b.i_fault_angle_deg,
                        "angle_shift_deg": float(metrics.angle_diff_deg(
                            b.i_fault_angle_deg, a.i_fault_angle_deg)),
                        "i2_over_i1_classical": a.negative_sequence_ratio,
                        "i2_over_i1_ibr": b.negative_sequence_ratio,
                        "iterations": b.iterations,
                        "converged": b.converged,
                        "well_posed": b.well_posed,
                        "note": b.note,
                    })
    return rows


def ieee2800_k2_sweep() -> list[dict]:
    """Sweep the negative-sequence injection gain K2. K2 = 0 is a pre-2800 inverter."""
    k2_values = config.scenarios()["ibr"]["ieee2800"]["k2_sweep"]
    rows = []
    for case in CASES:
        locations = config.scenarios()["faults"]["locations"].get(case, [])
        c, base, gens, seq = _setup(case, 60.0)
        if gens.size == 0:
            continue
        for bus in locations:
            for ft in ["SLG", "LL", "LLG"]:      # negative-sequence dependent
                a = F.solve_fault_classical(c, seq, bus, ft)
                for k2 in k2_values:
                    b = F.solve_fault_ibr_aware(c, seq, bus, ft, k2=k2)
                    rows.append({
                        "case": case, "bus": bus, "fault_type": ft, "k2": k2,
                        "if_pu": b.i_fault_mag,
                        "error_vs_classical_pct": metrics.pct_error(
                            b.i_fault_mag, a.i_fault_mag),
                        "i2_magnitude_pu": abs(b.i2),
                        "i2_over_i1": b.negative_sequence_ratio,
                        "converged": b.converged,
                        "well_posed": b.well_posed,
                    })
    return rows


def protection_study() -> list[dict]:
    """Relays set from a classical study, evaluated against IBR-aware faults.

    Relays are placed on every branch terminating at the faulted bus, which is
    where a real distance scheme would look first.
    """
    levels = [0, 40, 60, 80, 100]
    k2_values = [0.0, 4.0]                       # pre-2800 vs IEEE 2800
    rows = []

    for case in CASES:
        locations = config.scenarios()["faults"]["locations"].get(case, [])
        for pen in levels:
            c, base, gens, seq = _setup(case, pen)
            if gens.size == 0:
                continue
            br = c["branch"]
            for bus in locations:
                relays = [k for k in range(br.shape[0])
                          if int(br[k, 1]) == bus or int(br[k, 0]) == bus]
                for k in relays:
                    for ft in FAULT_TYPES:
                        a = F.solve_fault_classical(c, seq, bus, ft)
                        oa = PR.evaluate(c, seq, a, k, i_pickup_pu=2.0)
                        for k2 in k2_values:
                            b = F.solve_fault_ibr_aware(c, seq, bus, ft, k2=k2)
                            if not b.converged:
                                continue      # never score a relay on a bad solve
                            ob = PR.evaluate(c, seq, b, k, i_pickup_pu=2.0)
                            flags = PR.misoperation_flags(oa, ob)
                            rows.append({
                                "case": case,
                                "pen_pct": round(I.actual_penetration_pct(c, gens), 1),
                                "k2": k2,
                                "ieee2800": "yes" if k2 > 0 else "no",
                                "branch": f"{oa.from_bus}-{oa.to_bus}",
                                "fault_bus": bus,
                                "fault_type": ft,
                                "zone_classical": oa.zone,
                                "zone_ibr": ob.zone,
                                "direc_classical": oa.directional,
                                "direc_ibr": ob.directional,
                                "i2_classical": oa.i2_magnitude,
                                "i2_ibr": ob.i2_magnitude,
                                "t_oc_classical": oa.oc_operate_time_s,
                                "t_oc_ibr": ob.oc_operate_time_s,
                                "n_misoperations": sum(flags.values()),
                                **{f"flag_{k_}": v for k_, v in flags.items()},
                            })
    return rows


def main() -> None:
    config.ensure_output_dirs()
    print("=" * 78)
    print("WP3/WP4 -- fault analysis with current-limited IBRs, and protection")
    print("=" * 78)

    # ---------------------------------------------------------------- WP3a
    fe = fault_error_sweep()
    metrics.write_csv(config.TABLE_DIR / "wp3_fault_error.csv", fe)

    n_tot = len(fe)
    n_conv = sum(1 for r in fe if r["converged"])
    n_illposed = sum(1 for r in fe if not r["well_posed"])
    print(f"\nSolver status: {n_conv}/{n_tot} converged "
          f"({n_conv / n_tot * 100:.0f}%), {n_illposed} ill-posed.")
    print("  Non-convergence concentrates where synchronous sources run out. At")
    print("  100% IBR the phasor fault problem has no voltage reference at all and")
    print("  is reported as ill-posed rather than solved -- that boundary is a")
    print("  finding, not a defect. All statistics below use CONVERGED solves only.")

    print("\nFault current error, IBR-aware vs classical (converged solves only):")
    for case in CASES:
        pens = sorted({r["pen_pct"] for r in fe if r["case"] == case})
        print(f"\n  {case}   " + "".join(f"{ft:>13}" for ft in FAULT_TYPES))
        for pen in pens:
            cells = []
            for ft in FAULT_TYPES:
                vals = [r["error_pct"] for r in fe
                        if r["case"] == case and r["pen_pct"] == pen
                        and r["fault_type"] == ft and r["converged"]
                        and np.isfinite(r["error_pct"])]
                cells.append(f"{np.mean(vals):>12.2f}%" if vals else f"{'n/a':>13}")
            print(f"  {pen:5.1f}%" + "".join(cells))

    print("\nFault current ANGLE shift (deg) -- what directional elements see:")
    for case in CASES:
        pens = sorted({r["pen_pct"] for r in fe if r["case"] == case})
        print(f"\n  {case}   " + "".join(f"{ft:>13}" for ft in FAULT_TYPES))
        for pen in pens:
            cells = []
            for ft in FAULT_TYPES:
                vals = [r["angle_shift_deg"] for r in fe
                        if r["case"] == case and r["pen_pct"] == pen
                        and r["fault_type"] == ft and r["converged"]
                        and np.isfinite(r["angle_shift_deg"])]
                cells.append(f"{np.mean(vals):>+12.2f} " if vals else f"{'n/a':>13}")
            print(f"  {pen:5.1f}%" + "".join(cells))

    print("\nConvergence rate by penetration (phasor fault model validity):")
    for case in CASES:
        pens = sorted({r["pen_pct"] for r in fe if r["case"] == case})
        cells = []
        for pen in pens:
            sub = [r for r in fe if r["case"] == case and r["pen_pct"] == pen]
            cells.append(f"{sum(1 for r in sub if r['converged']) / len(sub) * 100:5.0f}%")
        print(f"  {case:>8}  " + "  ".join(f"{p:5.1f}%:{c}" for p, c in zip(pens, cells)))

    # ---------------------------------------------------------------- WP3b
    k2 = ieee2800_k2_sweep()
    metrics.write_csv(config.TABLE_DIR / "wp3_ieee2800_k2_sweep.csv", k2)
    print("\nIEEE 2800 negative-sequence injection: mean |I2| at the fault (pu)")
    print("  K2 = 0 is a pre-2800 inverter (no negative-sequence support)")
    ks = sorted({r["k2"] for r in k2})
    print(f"\n  {'case':>8} {'fault':>6} " + "".join(f"{f'K2={k:g}':>12}" for k in ks))
    for case in CASES:
        for ft in ["SLG", "LL", "LLG"]:
            cells = []
            for k in ks:
                vals = [r["i2_magnitude_pu"] for r in k2
                        if r["case"] == case and r["fault_type"] == ft
                        and r["k2"] == k and r["converged"]]
                cells.append(f"{np.mean(vals):>12.4f}" if vals else f"{'n/a':>12}")
            print(f"  {case:>8} {ft:>6} " + "".join(cells))

    # ----------------------------------------------------------------- WP4
    prot = protection_study()
    metrics.write_csv(config.TABLE_DIR / "wp4_protection_misoperation.csv", prot)

    print("\nProtection misoperations (relays set classically, faults solved IBR-aware):")
    print(f"  {'case':>8} {'pen%':>7} {'IEEE2800':>10} {'evaluated':>10} "
          f"{'misops':>8} {'rate':>8}")
    for case in CASES:
        pens = sorted({r["pen_pct"] for r in prot if r["case"] == case})
        for pen in pens:
            for tag in ["no", "yes"]:
                sub = [r for r in prot if r["case"] == case
                       and r["pen_pct"] == pen and r["ieee2800"] == tag]
                if not sub:
                    continue
                n_mis = sum(1 for r in sub if r["n_misoperations"] > 0)
                print(f"  {case:>8} {pen:7.1f} {tag:>10} {len(sub):>10} "
                      f"{n_mis:>8} {n_mis / len(sub) * 100:>7.1f}%")

    counts = Counter()
    for r in prot:
        for key, val in r.items():
            if key.startswith("flag_") and val:
                counts[key[5:]] += 1
    print("\n  Misoperation modes (all scenarios):")
    for name, n in counts.most_common():
        print(f"    {name:<24} {n:>5}")

    print(f"\nTables written to {config.TABLE_DIR}")


if __name__ == "__main__":
    main()
