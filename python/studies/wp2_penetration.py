"""WP2 -- IBR-aware power flow, system strength, and the penetration x loading plane.

The controlled experiment. Classical and IBR-aware power flow differ in exactly
one respect: where the reactive limit sits.

    classical:   Qmax = constant, read from the case file
    IBR-aware:   Qmax = sqrt((V * Ilim * Srating)^2 - P^2)

The inverter is sized so that at nominal voltage and rated power its capability
EQUALS the machine's case-file Qmax (verified to 1e-14 by `sizing_control`), so
the only thing that varies is how the limit MOVES with voltage and power.

Headline finding: at the nominal operating point the two agree exactly -- the
current limit does not bind. It binds under stress. So penetration is swept
AGAINST loading; sweeping it alone finds nothing, which is itself worth
reporting, because a study that swept penetration alone would conclude there is
no problem.

Outputs
    results/tables/wp2_sizing_control.csv
    results/tables/wp2_penetration_loading.csv
    results/tables/wp2_system_strength.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import (  # noqa: E402
    config, faults as F, ibr as I, metrics, ppc as P, solvers as S, strength as ST,
)

CASES = ["case14", "case30", "case39"]


def sizing_control() -> list[dict]:
    """Verify the experimental control: cap(V=1, P=PMAX) must equal case QMAX.

    If this does not hold, the comparison is measuring the sizing assumption
    rather than the physics, and every WP2 number is meaningless.
    """
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        for g in range(c["gen"].shape[0]):
            if c["is_slack_gen"][g] or c["gen"][g, 7] <= 0:
                continue
            s_rating = I.inverter_rating_pu(c, [g])
            pmax_pu = float(c["gen"][g, 8]) / c["baseMVA"]
            cap = I.q_capability_pu(1.0, pmax_pu, 1.0, s_rating) * c["baseMVA"]
            qmax = float(c["gen"][g, 3])
            rows.append({
                "case": case, "gen": g, "bus": int(c["gen"][g, 0]),
                "cap_mvar": cap, "case_qmax_mvar": qmax,
                "abs_error": abs(cap - qmax),
                "control_holds": bool(abs(cap - qmax) < 1e-6),
            })
    return rows


def penetration_loading_plane() -> list[dict]:
    """Sweep IBR penetration against load scaling; record where the models diverge."""
    levels = config.scenarios()["penetration"]["levels_pct"]
    loadings = config.scenarios()["ibr"]["load_scaling"]

    rows = []
    for case in CASES:
        c0 = P.load_ppc(case)
        for lam in loadings:
            c = I.scale_load(c0, lam)
            classical = S.newton_raphson(c, tol=1e-11, enforce_q_limits=True)
            for pen in levels:
                gens = I.select_ibr_gens(c, pen)
                r = I.ibr_powerflow(c, gens)

                row = {
                    "case": case,
                    "load_scale": lam,
                    "requested_pen_pct": pen,
                    "actual_pen_pct": r.penetration_pct,
                    "classical_converged": classical.converged,
                    "ibr_converged": r.converged,
                    "n_limited": len(r.limited_buses),
                    "limited_buses": " ".join(str(b) for b in r.limited_buses),
                }
                if classical.converged and r.converged:
                    row.update({
                        "max_dvm_pu": metrics.max_abs_error(r.vm, classical.vm),
                        "rmse_vm_pu": metrics.rmse(r.vm, classical.vm),
                        "vm_min_classical": float(classical.vm.min()),
                        "vm_min_ibr": float(r.vm.min()),
                    })
                else:
                    row.update({
                        "max_dvm_pu": float("nan"), "rmse_vm_pu": float("nan"),
                        "vm_min_classical": float("nan"), "vm_min_ibr": float("nan"),
                    })
                rows.append(row)
    return rows


def system_strength_sweep() -> list[dict]:
    """SCR / WSCR / CSCR at each penetration level."""
    levels = config.scenarios()["penetration"]["levels_pct"]
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        base = S.newton_raphson(c)
        for pen in levels:
            gens = I.select_ibr_gens(c, pen)
            seq = F.build_sequence_networks(
                c, base.v, ibr_buses=I.ibr_buses_from_gens(c, gens),
            )
            prof = ST.profile(c, seq, gens)
            rows.append({
                "case": case,
                "requested_pen_pct": pen,
                "actual_pen_pct": I.actual_penetration_pct(c, gens),
                "total_ibr_mw": prof["total_ibr_mw"],
                "min_scr": prof["min_scr"],
                "mean_scr": prof["mean_scr"],
                "wscr": prof["wscr"],
                "cscr": prof["cscr"],
                "n_weak_buses": prof["n_weak_buses"],
            })
    return rows


def main() -> None:
    config.ensure_output_dirs()
    print("=" * 78)
    print("WP2 -- IBR-aware power flow and system strength")
    print("=" * 78)

    ctrl = sizing_control()
    metrics.write_csv(config.TABLE_DIR / "wp2_sizing_control.csv", ctrl)
    bad = [r for r in ctrl if not r["control_holds"]]
    worst = max(r["abs_error"] for r in ctrl) if ctrl else float("nan")
    print(f"\nExperimental control -- inverter capability == machine QMAX at rated point:")
    print(f"  generators checked : {len(ctrl)}")
    print(f"  worst |error|      : {worst:.2e} MVAr")
    print(f"  control violations : {len(bad)}")
    if bad:
        print("  !! The comparison is NOT controlled; WP2 results are invalid.")

    plane = penetration_loading_plane()
    metrics.write_csv(config.TABLE_DIR / "wp2_penetration_loading.csv", plane)

    print("\nmax |dVm| between classical and IBR-aware power flow (pu):")
    levels = config.scenarios()["penetration"]["levels_pct"]
    for case in CASES:
        print(f"\n  {case}    " + "".join(f"{p:>10}%" for p in levels))
        for lam in config.scenarios()["ibr"]["load_scaling"]:
            cells = []
            for pen in levels:
                r = next(x for x in plane if x["case"] == case
                         and x["load_scale"] == lam and x["requested_pen_pct"] == pen)
                if not r["classical_converged"]:
                    cells.append(f"{'cl-div':>11}")
                elif not r["ibr_converged"]:
                    cells.append(f"{'ibr-div':>11}")
                else:
                    cells.append(f"{r['max_dvm_pu']:>11.5f}")
            print(f"  x{lam:<5.2f}" + "".join(cells))

    strength = system_strength_sweep()
    metrics.write_csv(config.TABLE_DIR / "wp2_system_strength.csv", strength)
    print("\nSystem strength (SCR >= 3.0 is the usual interconnection screen):")
    print(f"{'case':>8} {'pen%':>7} {'min SCR':>9} {'WSCR':>9} {'CSCR':>9} {'weak buses':>11}")
    for r in strength:
        print(f"{r['case']:>8} {r['actual_pen_pct']:7.1f} {r['min_scr']:9.2f} "
              f"{r['wscr']:9.2f} {r['cscr']:9.2f} {r['n_weak_buses']:11d}")

    print(f"\nTables written to {config.TABLE_DIR}")


if __name__ == "__main__":
    main()
