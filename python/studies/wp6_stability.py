"""WP6 -- Voltage stability, sensitivity, and N-1 security under IBR penetration.

Where WP2 asks "what is the operating point", this asks "how much margin is
left". That distinction matters more at high IBR penetration, not less:
displacing synchronous machines removes reactive reserve, and reactive reserve
is what voltage stability margin is made of.

Four studies:

  6a  CONTINUATION POWER FLOW -- loading margin (lambda_max) vs IBR penetration,
      with the converter capability curve applied to the displaced machines.

  6b  L-INDEX -- per-bus proximity to collapse from a single solved power flow,
      which is what an online monitoring scheme would actually use.

  6c  PTDF / LODF -- linear sensitivities, and a check that the fast DC-based
      N-1 screen agrees with the exact AC one.

  6d  N-1 CONTINGENCY SCREEN -- full AC branch-outage analysis, ranked, with
      islanding separated from genuine divergence.

Outputs
    results/tables/wp6_cpf_margin.csv
    results/tables/wp6_l_index.csv
    results/tables/wp6_sensitivity_check.csv
    results/tables/wp6_n1_ranking.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import (  # noqa: E402
    config, ibr as I, metrics, ppc as P, solvers as S, stability as ST, ybus as Y,
)

CASES = ["case14", "case30", "case39"]


def _apply_converter_limits(ppc: dict, ibr_gens: np.ndarray) -> dict:
    """Replace the case-file Qmax on displaced machines with converter capability.

    Same substitution WP2 makes, so the CPF is measuring the same model change
    as the power-flow study rather than a different one.
    """
    out = {**ppc, "gen": ppc["gen"].copy()}
    for g in np.asarray(ibr_gens, dtype=int):
        s_rating = I.inverter_rating_pu(ppc, [int(g)])
        p_pu = float(ppc["gen"][g, P.PG]) / ppc["baseMVA"]
        cap = I.q_capability_pu(1.0, p_pu, 1.0, s_rating) * ppc["baseMVA"]
        out["gen"][g, P.QMAX] = cap
        out["gen"][g, P.QMIN] = -cap
    return out


def cpf_margin_sweep() -> list[dict]:
    """Loading margin vs IBR penetration."""
    levels = config.scenarios()["penetration"]["levels_pct"]
    rows = []
    for case in CASES:
        c0 = P.load_ppc(case)
        for pen in levels:
            gens = I.select_ibr_gens(c0, pen)
            c = _apply_converter_limits(c0, gens) if gens.size else c0
            r = ST.continuation_power_flow(c)
            rows.append({
                "case": case,
                "requested_pen_pct": pen,
                "actual_pen_pct": I.actual_penetration_pct(c0, gens),
                "lambda_max": r.lambda_max,
                "loading_margin_pct": r.loading_margin_pct,
                "critical_bus": r.critical_bus,
                "n_points": len(r.lambdas),
                "converged": r.converged,
            })
    return rows


def l_index_sweep() -> list[dict]:
    """Worst-bus L-index vs penetration. 0 = no load, 1 = collapse."""
    levels = config.scenarios()["penetration"]["levels_pct"]
    loadings = [1.0, 1.2, 1.4]
    rows = []
    for case in CASES:
        c0 = P.load_ppc(case)
        for lam in loadings:
            base = I.scale_load(c0, lam)
            for pen in levels:
                gens = I.select_ibr_gens(base, pen)
                c = _apply_converter_limits(base, gens) if gens.size else base
                pf = S.newton_raphson(c, tol=1e-10, enforce_q_limits=True)
                if not pf.converged:
                    continue
                li = ST.l_index(c, pf.v)
                finite = li[np.isfinite(li)]
                if finite.size == 0:
                    continue
                worst = int(np.nanargmax(li))
                rows.append({
                    "case": case,
                    "load_scale": lam,
                    "pen_pct": I.actual_penetration_pct(base, gens),
                    "max_l_index": float(np.nanmax(li)),
                    "mean_l_index": float(np.mean(finite)),
                    "worst_bus": int(c["bus"][worst, P.BUS_I]),
                    "vm_min_pu": float(pf.vm.min()),
                })
    return rows


def sensitivity_check() -> list[dict]:
    """Validate the linear (DC) N-1 screen against the exact AC one.

    LODF gives every branch outage from one matrix, which is how real-time
    contingency analysis keeps up with the clock. It is only useful if it agrees
    with the truth, so this measures the disagreement rather than assuming it.
    """
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        pf = S.newton_raphson(c)
        h = ST.ptdf(c)
        ld = ST.lodf(c, h)

        sf, _ = Y.branch_flows(c, pf.v)
        base_flow = np.real(sf) * c["baseMVA"]

        ac = {r.branch: r for r in ST.n1_screen(c)}
        for k, res in ac.items():
            if not res.converged:
                continue
            # DC prediction of post-outage flows
            predicted = base_flow + ld[:, k] * base_flow[k]
            predicted[k] = 0.0

            trial = {**c, "branch": c["branch"].copy()}
            trial["branch"][k, P.BR_STATUS] = 0
            pf_out = S.newton_raphson(trial, tol=1e-9, enforce_q_limits=True)
            if not pf_out.converged:
                continue
            sf_out, _ = Y.branch_flows(trial, pf_out.v)
            actual = np.real(sf_out) * c["baseMVA"]
            actual[k] = 0.0

            rows.append({
                "case": case,
                "outaged_branch": f"{int(c['branch'][k, P.F_BUS])}-{int(c['branch'][k, P.T_BUS])}",
                "max_flow_error_mw": metrics.max_abs_error(predicted, actual),
                "rmse_flow_mw": metrics.rmse(predicted, actual),
                "max_actual_flow_mw": float(np.max(np.abs(actual))),
            })
    return rows


def n1_ranking() -> list[dict]:
    """Full AC N-1 screen at 0% and high IBR penetration."""
    rows = []
    for case in CASES:
        c0 = P.load_ppc(case)
        for pen in [0, 80]:
            gens = I.select_ibr_gens(c0, pen)
            c = _apply_converter_limits(c0, gens) if gens.size else c0
            for r in ST.n1_screen(c):
                rows.append({
                    "case": case,
                    "pen_pct": round(I.actual_penetration_pct(c0, gens), 1),
                    "branch": f"{r.from_bus}-{r.to_bus}",
                    "status": r.status,
                    "n_violations": r.n_violations,
                    "max_loading_pct": r.max_loading_pct,
                    "min_vm_pu": r.min_vm_pu,
                })
    return rows


def main() -> None:
    config.ensure_output_dirs()
    print("=" * 78)
    print("WP6 -- voltage stability, sensitivity and N-1 security")
    print("=" * 78)

    cpf = cpf_margin_sweep()
    metrics.write_csv(config.TABLE_DIR / "wp6_cpf_margin.csv", cpf)
    print("\nLoading margin from continuation power flow (% load increase to collapse):")
    levels = config.scenarios()["penetration"]["levels_pct"]
    print(f"  {'case':>8} " + "".join(f"{p:>10}%" for p in levels))
    for case in CASES:
        cells = []
        for pen in levels:
            r = next((x for x in cpf if x["case"] == case
                      and x["requested_pen_pct"] == pen), None)
            cells.append(f"{r['loading_margin_pct']:10.1f} " if r else f"{'n/a':>11}")
        print(f"  {case:>8} " + "".join(cells))
    print("  (critical bus: " + ", ".join(
        f"{c}={next(x['critical_bus'] for x in cpf if x['case'] == c)}" for c in CASES) + ")")

    li = l_index_sweep()
    metrics.write_csv(config.TABLE_DIR / "wp6_l_index.csv", li)
    print("\nWorst-bus L-index (0 = unloaded, 1 = voltage collapse):")
    print(f"  {'case':>8} {'load':>6} " + "".join(f"{p:>9}%" for p in levels))
    for case in CASES:
        for lam in [1.0, 1.2, 1.4]:
            cells = []
            for pen in levels:
                r = next((x for x in li if x["case"] == case and x["load_scale"] == lam
                          and abs(x["pen_pct"] - I.actual_penetration_pct(
                              P.load_ppc(case), I.select_ibr_gens(P.load_ppc(case), pen))) < 0.1), None)
                cells.append(f"{r['max_l_index']:10.4f}" if r else f"{'---':>10}")
            print(f"  {case:>8} x{lam:<5.2f}" + "".join(cells))

    sc = sensitivity_check()
    metrics.write_csv(config.TABLE_DIR / "wp6_sensitivity_check.csv", sc)
    print("\nLinear (LODF) N-1 screen vs exact AC, post-outage branch flows:")
    print(f"  {'case':>8} {'outages':>8} {'max error':>12} {'mean RMSE':>12} {'as % of flow':>13}")
    for case in CASES:
        sub = [r for r in sc if r["case"] == case]
        if not sub:
            continue
        worst = max(r["max_flow_error_mw"] for r in sub)
        mean_rmse = float(np.mean([r["rmse_flow_mw"] for r in sub]))
        peak = max(r["max_actual_flow_mw"] for r in sub)
        print(f"  {case:>8} {len(sub):>8} {worst:>11.2f}M {mean_rmse:>11.2f}M "
              f"{worst / peak * 100:>12.1f}%")

    n1 = n1_ranking()
    metrics.write_csv(config.TABLE_DIR / "wp6_n1_ranking.csv", n1)
    print("\nN-1 branch-outage screen:")
    print(f"  {'case':>8} {'pen%':>6} {'solved':>7} {'islands':>8} {'diverged':>9} "
          f"{'with violations':>16}")
    for case in CASES:
        for pen in sorted({r["pen_pct"] for r in n1 if r["case"] == case}):
            sub = [r for r in n1 if r["case"] == case and r["pen_pct"] == pen]
            solved = [r for r in sub if r["status"] == "solved"]
            print(f"  {case:>8} {pen:6.1f} {len(solved):>7} "
                  f"{sum(1 for r in sub if r['status'] == 'islands network'):>8} "
                  f"{sum(1 for r in sub if r['status'] == 'diverged'):>9} "
                  f"{sum(1 for r in solved if r['n_violations'] > 0):>16}")

    print("\n  Worst SOLVED contingencies by new violations (IEEE 39-bus, 0% IBR):")
    solved = [r for r in n1 if r["case"] == "case39" and r["pen_pct"] == 0.0
              and r["status"] == "solved"]
    solved.sort(key=lambda r: (-int(float(r["n_violations"])),
                               -(float(r["max_loading_pct"])
                                 if r["max_loading_pct"] not in ("", "nan") else 0.0)))
    for r in solved[:5]:
        print(f"    {r['branch']:>8}  new violations={r['n_violations']:<3} "
              f"max loading={float(r['max_loading_pct']):6.1f}%  "
              f"Vmin={float(r['min_vm_pu']):.4f}")

    print(f"\nTables written to {config.TABLE_DIR}")


if __name__ == "__main__":
    main()
