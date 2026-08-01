"""WP7 -- Transient stability, critical clearing time, and a metric that breaks.

Three parts, and the middle one is a negative result that matters:

  7a  BASELINE VALIDATION. CCT on the intact system, checked against published
      values for these benchmarks. If this is wrong nothing else here counts.

  7b  THE METRIC BREAKS. Naively sweeping grid-following penetration makes CCT
      look BETTER, which is false. Converting a machine to a grid-following
      converter removes its rotor from the swing model, so the rotor-angle
      criterion loses the failure mode it exists to detect. Reported explicitly,
      because this exact mistake is easy to make and would invert the
      conclusion of the whole project.

  7c  THE CONTROLLED TEST. Holding topology, dispatch and unit count fixed and
      varying ONLY grid-forming virtual inertia. This is the comparison the
      model can answer honestly, and it is where `gfm_extends_cct` is settled.

Outputs
    results/tables/wp7_cct_baseline.csv
    results/tables/wp7_cct_metric_artefact.csv
    results/tables/wp7_virtual_inertia_sweep.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import (  # noqa: E402
    config, ibr as I, metrics, ppc as P, solvers as S, transient as T,
)

# fault bus -> the adjacent line tripped to clear it
SCENARIOS = {
    "case9": (7, (6, 7)),
    "case39": (16, (16, 17)),
}
PUBLISHED_CCT = {"case39": "~0.15 s (Athay et al., 3LG at bus 16)"}


def _branch_index(ppc, pair) -> int | None:
    br = ppc["branch"]
    for k in range(br.shape[0]):
        if {int(br[k, 0]), int(br[k, 1])} == set(pair):
            return k
    return None


def baseline_cct() -> list[dict]:
    rows = []
    for case, (fbus, pair) in SCENARIOS.items():
        c = P.load_ppc(case)
        pf = S.newton_raphson(c)
        k = _branch_index(c, pair)
        net = T.build_dynamic_network(c, pf.v)
        net.damping[:] = 0.0
        r = T.critical_clearing_time(c, pf.v, net, fbus, outaged_branch=k, t_max=1.0)
        rows.append({
            "case": case,
            "fault_bus": fbus,
            "cleared_by_tripping": f"{pair[0]}-{pair[1]}",
            "n_machines": net.n,
            "h_total_s": float(net.h.sum()),
            "cct_s": r["cct_s"],
            "bracketed": r["bracketed"],
            "published_reference": PUBLISHED_CCT.get(case, "n/a"),
        })
    return rows


def metric_artefact() -> list[dict]:
    """Demonstrate that CCT is NOT comparable across differing unit counts."""
    rows = []
    case = "case39"
    fbus, pair = SCENARIOS[case]
    c = P.load_ppc(case)
    pf = S.newton_raphson(c)
    k = _branch_index(c, pair)

    for pen in [0, 40, 60, 80]:
        gens = I.select_ibr_gens(c, pen)
        net = T.build_dynamic_network(c, pf.v, ibr_gens=gens)
        net.damping[:] = 0.0
        r = T.critical_clearing_time(
            c, pf.v, net, fbus, ibr_gens=gens, outaged_branch=k, t_max=1.0)
        rows.append({
            "case": case,
            "pen_pct": round(I.actual_penetration_pct(c, gens), 1),
            "gfm_share_pct": 0,
            "n_dynamic_units": net.n,
            "h_total_s": float(net.h.sum()),
            "cct_s": r["cct_s"],
            "valid": pen == 0,
            "why": "" if pen == 0 else
                   "NOT COMPARABLE: fewer rotors, so the rotor-angle criterion "
                   "no longer measures the same failure mode",
        })
    return rows


def virtual_inertia() -> list[dict]:
    """The controlled comparison: only virtual inertia varies."""
    h_values = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    rows = []
    for case, (fbus, pair) in SCENARIOS.items():
        c = P.load_ppc(case)
        pf = S.newton_raphson(c)
        k = _branch_index(c, pair)
        seen: set[float] = set()
        for pen in [60, 80]:
            gens = I.select_ibr_gens(c, pen)
            if gens.size == 0:
                continue
            # small cases hit their penetration ceiling, so two requested levels
            # can resolve to the same fleet -- run it once
            actual = round(I.actual_penetration_pct(c, gens), 1)
            if actual in seen:
                continue
            seen.add(actual)
            _, gfm = I.split_gfl_gfm(c, gens, 100.0)
            for r in T.virtual_inertia_sweep(
                    c, pf.v, fbus, gens, gfm, h_values, outaged_branch=k, t_max=1.0):
                rows.append({
                    "case": case,
                    "pen_pct": round(I.actual_penetration_pct(c, gens), 1),
                    "gfm_share_pct": 100,
                    **r,
                })
    return rows


def main() -> None:
    config.ensure_output_dirs()
    print("=" * 78)
    print("WP7 -- transient stability and critical clearing time")
    print("=" * 78)

    base = baseline_cct()
    metrics.write_csv(config.TABLE_DIR / "wp7_cct_baseline.csv", base)
    print("\n7a. Baseline CCT on the intact system (all synchronous):")
    print(f"  {'case':>8} {'fault':>6} {'cleared by':>12} {'units':>6} "
          f"{'H_total':>9} {'CCT (s)':>9}   reference")
    for r in base:
        print(f"  {r['case']:>8} {r['fault_bus']:>6} {r['cleared_by_tripping']:>12} "
              f"{r['n_machines']:>6} {r['h_total_s']:>9.1f} {r['cct_s']:>9.3f}   "
              f"{r['published_reference']}")

    art = metric_artefact()
    metrics.write_csv(config.TABLE_DIR / "wp7_cct_metric_artefact.csv", art)
    print("\n7b. Why a naive GFL penetration sweep is INVALID:")
    print(f"  {'pen%':>7} {'rotors':>8} {'H_total':>9} {'CCT (s)':>9}  verdict")
    for r in art:
        verdict = "valid baseline" if r["valid"] else "ARTEFACT -- do not report"
        print(f"  {r['pen_pct']:>7.1f} {r['n_dynamic_units']:>8} "
              f"{r['h_total_s']:>9.1f} {r['cct_s']:>9.3f}  {verdict}")
    print("\n  CCT appears to IMPROVE as grid-following penetration rises. It has not.")
    print("  Converting a machine to a grid-following converter deletes its rotor,")
    print("  and the stability criterion is rotor-angle separation -- so the metric")
    print("  loses the failure mode it exists to detect. A grid-following converter")
    print("  destabilises through PLL loss of synchronisation, which is a converter")
    print("  control phenomenon an electromechanical model cannot represent at all.")
    print("  This is an independent demonstration that EMT analysis is mandatory.")

    vi = virtual_inertia()
    metrics.write_csv(config.TABLE_DIR / "wp7_virtual_inertia_sweep.csv", vi)
    print("\n7c. CONTROLLED test -- topology, dispatch and unit count all fixed,")
    print("    only grid-forming virtual inertia varies:")
    for case in sorted({r["case"] for r in vi}):
        for pen in sorted({r["pen_pct"] for r in vi if r["case"] == case}):
            sub = [r for r in vi if r["case"] == case and r["pen_pct"] == pen]
            sub.sort(key=lambda r: r["h_virtual_s"])
            if not sub:
                continue
            print(f"\n    {case} @ {pen}% IBR, 100% grid-forming "
                  f"({sub[0]['n_dynamic_units']} dynamic units throughout)")
            print(f"      {'H_virtual (s)':>14} {'H_total':>9} {'CCT (s)':>9}")
            for r in sub:
                print(f"      {r['h_virtual_s']:>14.1f} {r['h_total_s']:>9.1f} "
                      f"{r['cct_s']:>9.3f}")
            lo, hi = sub[0]["cct_s"], sub[-1]["cct_s"]
            if lo > 0:
                print(f"      -> CCT {lo:.3f}s -> {hi:.3f}s "
                      f"({hi / lo:.1f}x) across the inertia range")

    print(f"\nTables written to {config.TABLE_DIR}")


if __name__ == "__main__":
    main()
