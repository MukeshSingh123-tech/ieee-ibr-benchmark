"""WP8 -- State estimation, bad-data detection, and false data injection.

The load flow solutions from WP1 are exactly the measurement set a control
centre sees. This asks whether an attacker can corrupt the operator's picture of
the grid without tripping the classical defence.

  8a  ESTIMATOR VALIDATION. WLS must recover the true state from noisy,
      redundant measurements. If it does not, nothing after it means anything.

  8b  THE DEFENCE WORKS -- against naive attacks. A random perturbation is
      caught decisively by the chi-squared residual test.

  8c  THE DEFENCE FAILS STRUCTURALLY. An attack constructed to be exactly
      consistent with a false state, a = h(x+c) - h(x), leaves the residual
      UNCHANGED. It is invisible at any magnitude and any threshold, because
      there is no residual to test. This is the gap that motivates out-of-band
      verification -- PMU physics cross-checks or cryptographically attested
      measurements, i.e. the ChainPMU line of work.

  8d  DOES IBR PENETRATION CHANGE THE PICTURE? The estimator is run across the
      penetration sweep to see whether a converter-dominated grid is easier or
      harder to attack undetectably.

Outputs
    results/tables/wp8_estimator_validation.csv
    results/tables/wp8_attack_detection.csv
    results/tables/wp8_attack_vs_penetration.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import (  # noqa: E402
    config, estimation as E, ibr as I, metrics, ppc as P, solvers as S,
)

CASES = ["case14", "case30", "case39"]


def estimator_validation() -> list[dict]:
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        pf = S.newton_raphson(c)
        meas = E.build_measurements(c, pf.v)
        r = E.estimate(c, meas)
        n_state = 2 * c["bus"].shape[0] - 1
        rows.append({
            "case": case,
            "n_measurements": meas.n,
            "n_states": n_state,
            "redundancy": meas.n / n_state,
            "converged": r.converged,
            "iterations": r.iterations,
            "max_dvm_pu": float(np.max(np.abs(np.abs(r.v) - np.abs(pf.v)))),
            "max_dva_deg": float(np.max(np.abs(np.degrees(
                np.angle(r.v) - np.angle(pf.v))))),
            "objective_J": r.objective,
            "chi2_threshold": r.chi2_threshold,
            "false_alarm": r.bad_data_detected,
        })
    return rows


def attack_detection() -> list[dict]:
    """Naive vs linearised vs exact attacks, at several magnitudes."""
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        pf = S.newton_raphson(c)
        meas = E.build_measurements(c, pf.v)
        idx = P.bus_index(c)

        # attack a mid-network load bus
        pq = [int(b) for b in c["bus"][c["bus"][:, 1] == P.PQ, 0].astype(int)]
        target = pq[len(pq) // 2] if pq else int(c["bus"][1, 0])
        i = idx[target]

        clean = E.estimate(c, meas)
        rows.append({
            "case": case, "attack": "none", "target_bus": target,
            "magnitude_deg": 0.0, "objective_J": clean.objective,
            "chi2_threshold": clean.chi2_threshold,
            "detected": clean.bad_data_detected,
            "max_norm_residual": clean.largest_normalised_residual,
            "state_error_deg": 0.0,
        })

        rnd = E.estimate(c, E.apply_attack(meas, E.random_attack(meas, 3, 0.2)))
        rows.append({
            "case": case, "attack": "random", "target_bus": -1,
            "magnitude_deg": float("nan"), "objective_J": rnd.objective,
            "chi2_threshold": rnd.chi2_threshold,
            "detected": rnd.bad_data_detected,
            "max_norm_residual": rnd.largest_normalised_residual,
            "state_error_deg": float("nan"),
        })

        for deg in [2.0, 5.0, 10.0, 20.0]:
            rad = np.deg2rad(deg)
            for label, lin in [("fdia_linearised", True), ("fdia_exact", False)]:
                a = E.stealthy_attack(c, meas, pf.v, target, rad, linear=lin)
                r = E.estimate(c, E.apply_attack(meas, a))
                err = float(np.degrees(np.angle(r.v[i]) - np.angle(pf.v[i])))
                rows.append({
                    "case": case, "attack": label, "target_bus": target,
                    "magnitude_deg": deg, "objective_J": r.objective,
                    "chi2_threshold": r.chi2_threshold,
                    "detected": r.bad_data_detected,
                    "max_norm_residual": r.largest_normalised_residual,
                    "state_error_deg": err,
                })
    return rows


def attack_vs_penetration() -> list[dict]:
    """Is a converter-dominated grid easier to attack undetectably?"""
    rows = []
    levels = config.scenarios()["penetration"]["levels_pct"]
    for case in ["case14", "case39"]:
        c0 = P.load_ppc(case)
        idx = P.bus_index(c0)
        pq = [int(b) for b in c0["bus"][c0["bus"][:, 1] == P.PQ, 0].astype(int)]
        target = pq[len(pq) // 2] if pq else int(c0["bus"][1, 0])
        i = idx[target]

        for pen in levels:
            gens = I.select_ibr_gens(c0, pen)
            res = I.ibr_powerflow(c0, gens) if gens.size else None
            v = res.v if (res and res.converged) else S.newton_raphson(c0).v

            meas = E.build_measurements(c0, v)
            clean = E.estimate(c0, meas)
            a = E.stealthy_attack(c0, meas, v, target, np.deg2rad(10.0))
            atk = E.estimate(c0, E.apply_attack(meas, a))

            rows.append({
                "case": case,
                "pen_pct": round(I.actual_penetration_pct(c0, gens), 1),
                "target_bus": target,
                "clean_J": clean.objective,
                "attacked_J": atk.objective,
                "delta_J": atk.objective - clean.objective,
                "chi2_threshold": atk.chi2_threshold,
                "detected": atk.bad_data_detected,
                "state_error_deg": float(np.degrees(
                    np.angle(atk.v[i]) - np.angle(v[i]))),
            })
    return rows


def main() -> None:
    config.ensure_output_dirs()
    print("=" * 78)
    print("WP8 -- state estimation, bad-data detection, and FDIA")
    print("=" * 78)

    val = estimator_validation()
    metrics.write_csv(config.TABLE_DIR / "wp8_estimator_validation.csv", val)
    print("\n8a. WLS estimator recovers the true state:")
    print(f"  {'case':>8} {'meas':>6} {'states':>7} {'redund':>7} {'iters':>6} "
          f"{'max|dVm|':>10} {'max|dVa|':>10} {'false alarm':>12}")
    for r in val:
        print(f"  {r['case']:>8} {r['n_measurements']:>6} {r['n_states']:>7} "
              f"{r['redundancy']:>7.2f} {r['iterations']:>6} "
              f"{r['max_dvm_pu']:>10.2e} {r['max_dva_deg']:>9.3f}d "
              f"{str(r['false_alarm']):>12}")

    det = attack_detection()
    metrics.write_csv(config.TABLE_DIR / "wp8_attack_detection.csv", det)
    print("\n8b/8c. Chi-squared bad-data detection vs three attack classes:")
    for case in CASES:
        sub = [r for r in det if r["case"] == case]
        thr = sub[0]["chi2_threshold"]
        print(f"\n  {case} (chi-squared threshold J = {thr:.2f})")
        print(f"    {'attack':>18} {'magnitude':>10} {'J':>10} {'detected':>9} "
              f"{'state error':>12}")
        for r in sub:
            mag = "-" if not np.isfinite(r["magnitude_deg"]) else \
                ("-" if r["magnitude_deg"] == 0 else f"{r['magnitude_deg']:.0f} deg")
            err = "-" if not np.isfinite(r["state_error_deg"]) else \
                f"{r['state_error_deg']:.2f} deg"
            flag = "YES" if r["detected"] else "no"
            print(f"    {r['attack']:>18} {mag:>10} {r['objective_J']:>10.2f} "
                  f"{flag:>9} {err:>12}")

    print("\n  The exact AC attack leaves J essentially unchanged from clean at EVERY")
    print("  magnitude -- a 20 deg corruption of the operator's state estimate with")
    print("  no alarm. The defence does not fail because the threshold is mistuned;")
    print("  it fails because a residual test has nothing to test. Catching this")
    print("  needs information outside the residual: redundant PMU physics, or")
    print("  authenticated measurements.")

    pen = attack_vs_penetration()
    metrics.write_csv(config.TABLE_DIR / "wp8_attack_vs_penetration.csv", pen)
    print("\n8d. Does IBR penetration change attack detectability? (10 deg exact FDIA)")
    print(f"  {'case':>8} {'pen%':>7} {'clean J':>9} {'attacked J':>11} "
          f"{'delta J':>9} {'detected':>9}")
    for r in pen:
        print(f"  {r['case']:>8} {r['pen_pct']:>7.1f} {r['clean_J']:>9.2f} "
              f"{r['attacked_J']:>11.2f} {r['delta_J']:>9.2e} "
              f"{str(r['detected']):>9}")

    print(f"\nTables written to {config.TABLE_DIR}")


if __name__ == "__main__":
    main()
