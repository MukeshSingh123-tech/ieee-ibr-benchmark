"""WP1 -- Classical load flow baseline and solver comparison.

This is the original project brief ("Newton-Raphson / Gauss-Seidel load flow ...
with a comparison of results and convergence"), executed properly and used as
the CONTROL GROUP for everything that follows. Two things must come out of it:

  1. Every solver agrees on the answer, and agrees with pandapower. If they do
     not, no later result can be trusted.
  2. They disagree wildly on the COST of getting there, and that disagreement
     scales with system size -- which is the actual content of the comparison.

Outputs
    results/tables/wp1_convergence.csv
    results/tables/wp1_crosstool.csv
    results/tables/wp1_gs_acceleration.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import config, metrics, ppc as P, solvers as S, ybus as Y  # noqa: E402

CASES = ["case9", "case14", "case30", "case39"]
ALGORITHMS = ["nr", "gs", "fdxb", "fdbx"]


def convergence_comparison() -> list[dict]:
    """Iterations, time, and final mismatch for each solver on each case."""
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        for alg in ALGORITHMS:
            r = S.solve(c, alg)
            rows.append({
                "case": case,
                "n_bus": c["bus"].shape[0],
                "algorithm": alg,
                "converged": r.converged,
                "iterations": r.iterations,
                "time_ms": r.elapsed_s * 1e3,
                "final_mismatch": r.final_mismatch,
                "jacobian_cond": r.jacobian_cond,
            })
    return rows


def cross_tool_validation() -> list[dict]:
    """Hand-written solvers vs pandapower -- the gate that licenses everything else."""
    import pandapower as pp
    import pandapower.networks as pn

    def gate_for(alg: str) -> tuple[float, float]:
        key = "handwritten_gs_vs_matpower" if alg == "gs" else "handwritten_vs_matpower"
        return (config.tol("cross_tool", key, "vm_pu"),
                config.tol("cross_tool", key, "va_deg"))

    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        net = getattr(pn, case)()
        pp.runpp(net, numba=False)
        ref_vm = net.res_bus.vm_pu.values
        ref_va = net.res_bus.va_degree.values

        ybus_err = Y.validate_against_reference(c)

        for alg in ALGORITHMS:
            tol_vm, tol_va = gate_for(alg)
            r = S.solve(c, alg)
            if not r.converged:
                rows.append({
                    "case": case, "algorithm": alg, "converged": False,
                    "max_dvm_pu": float("nan"), "max_dva_deg": float("nan"),
                    "ybus_err": ybus_err, "gate_vm_pu": tol_vm, "passes_gate": False,
                })
                continue
            d_vm = metrics.max_abs_error(r.vm, ref_vm)
            d_va = float(np.max(np.abs(metrics.angle_diff_deg(r.va_deg, ref_va))))
            rows.append({
                "case": case, "algorithm": alg, "converged": True,
                "max_dvm_pu": d_vm, "max_dva_deg": d_va,
                "ybus_err": ybus_err, "gate_vm_pu": tol_vm,
                "passes_gate": bool(d_vm <= tol_vm and d_va <= tol_va),
            })
    return rows


def gauss_seidel_acceleration() -> list[dict]:
    """Characterise Gauss-Seidel rather than just recording that it failed.

    Plain GS on IEEE 39-bus does not converge within 5000 iterations at the
    default acceleration factor. Reporting "diverged" would be true but useless;
    sweeping the acceleration factor shows WHY -- GS convergence depends on a
    parameter that Newton-Raphson does not have and does not need.
    """
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        for accel in [1.0, 1.2, 1.4, 1.6, 1.8]:
            r = S.gauss_seidel(c, accel=accel, max_iter=8000)
            rows.append({
                "case": case,
                "accel": accel,
                "converged": r.converged,
                "iterations": r.iterations,
                "time_ms": r.elapsed_s * 1e3,
                "final_mismatch": r.final_mismatch,
            })
    return rows


def main() -> None:
    config.ensure_output_dirs()

    print("=" * 78)
    print("WP1 -- classical load flow baseline")
    print("=" * 78)

    conv = convergence_comparison()
    metrics.write_csv(config.TABLE_DIR / "wp1_convergence.csv", conv)
    print("\nSolver convergence (iterations to 1e-10 mismatch):")
    print(f"{'case':>8} {'n_bus':>6} " + "".join(f"{a:>12}" for a in ALGORITHMS))
    for case in CASES:
        cells = []
        for alg in ALGORITHMS:
            r = next(x for x in conv if x["case"] == case and x["algorithm"] == alg)
            cells.append(f"{r['iterations']:>12}" if r["converged"] else f"{'DIVERGED':>12}")
        n_bus = next(x["n_bus"] for x in conv if x["case"] == case)
        print(f"{case:>8} {n_bus:>6} " + "".join(cells))

    cross = cross_tool_validation()
    metrics.write_csv(config.TABLE_DIR / "wp1_crosstool.csv", cross)
    failures = [r for r in cross if not r["passes_gate"]]
    print(f"\nCross-tool gate vs pandapower "
          f"(tol {config.tol('cross_tool', 'handwritten_vs_matpower', 'vm_pu'):.0e} pu):")
    worst = max((r["max_dvm_pu"] for r in cross if np.isfinite(r["max_dvm_pu"])), default=float("nan"))
    print(f"  worst |dVm| across all solvers/cases : {worst:.3e} pu")
    print(f"  Ybus vs pandapower internal          : "
          f"{max(r['ybus_err'] for r in cross):.3e}")
    print(f"  gate failures                        : {len(failures)}")
    for r in failures:
        print(f"    FAIL {r['case']}/{r['algorithm']}: "
              f"dVm={r['max_dvm_pu']:.3e} vs gate {r['gate_vm_pu']:.0e}"
              + ("  (solver did not converge)" if not r["converged"] else ""))

    accel = gauss_seidel_acceleration()
    metrics.write_csv(config.TABLE_DIR / "wp1_gs_acceleration.csv", accel)
    print("\nGauss-Seidel iterations vs acceleration factor:")
    print(f"{'case':>8} " + "".join(f"{a:>10}" for a in [1.0, 1.2, 1.4, 1.6, 1.8]))
    for case in CASES:
        cells = []
        for a in [1.0, 1.2, 1.4, 1.6, 1.8]:
            r = next(x for x in accel if x["case"] == case and x["accel"] == a)
            cells.append(f"{r['iterations']:>10}" if r["converged"] else f"{'---':>10}")
        print(f"{case:>8} " + "".join(cells))

    print(f"\nTables written to {config.TABLE_DIR}")


if __name__ == "__main__":
    main()
