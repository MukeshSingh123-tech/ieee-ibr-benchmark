"""Generate every figure in the report from the CSV tables in results/tables/.

Deliberately reads the CSVs rather than recomputing: a figure that cannot be
regenerated from a committed table is a figure whose provenance is unclear, and
this way the numbers in the report and the numbers in the figures cannot drift.

Run the studies first, then:
    python studies/make_figures.py
"""

from __future__ import annotations

import csv
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import config, plotting as V  # noqa: E402


def load(name: str) -> list[dict]:
    path = config.TABLE_DIR / name
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def fnum(row: dict, key: str) -> float:
    try:
        v = row.get(key, "")
        return float(v) if v not in ("", "nan", "None") else float("nan")
    except (TypeError, ValueError):
        return float("nan")


# =============================================================================
def fig_wp1_convergence() -> None:
    rows = load("wp1_convergence.csv")
    if not rows:
        return
    cases = ["case9", "case14", "case30", "case39"]
    algs = ["nr", "fdxb", "fdbx", "gs"]
    names = {"nr": "Newton-Raphson", "fdxb": "Fast-decoupled XB",
             "fdbx": "Fast-decoupled BX", "gs": "Gauss-Seidel"}

    series = {}
    for a in algs:
        vals = []
        for c in cases:
            r = next((x for x in rows if x["case"] == c and x["algorithm"] == a), None)
            vals.append(fnum(r, "iterations") if r and r["converged"] == "True"
                        else float("nan"))
        series[names[a]] = vals

    fig = V.grouped_bars(
        [c.replace("case", "IEEE ") + "-bus" for c in cases], series,
        title="Newton-Raphson cost is flat in system size; Gauss-Seidel is not",
        xlabel="Test system", ylabel="Iterations to 1e-10 mismatch",
        logy=True,
        note="Gauss-Seidel bar absent for IEEE 39-bus: it does not converge within "
             "5000 iterations at acceleration 1.6.")
    print("  ", V.save(fig, "wp1_convergence"))


def fig_wp1_gs_acceleration() -> None:
    rows = load("wp1_gs_acceleration.csv")
    if not rows:
        return
    cases = ["case9", "case14", "case30", "case39"]
    accels = sorted({fnum(r, "accel") for r in rows})

    series = {}
    for c in cases:
        vals = []
        for a in accels:
            r = next((x for x in rows if x["case"] == c
                      and abs(fnum(x, "accel") - a) < 1e-9), None)
            vals.append(fnum(r, "iterations") if r and r["converged"] == "True"
                        else float("nan"))
        series[c.replace("case", "IEEE ") + "-bus"] = vals

    fig = V.line_by_series(
        accels, series,
        title="Gauss-Seidel depends on a tuning parameter Newton-Raphson does not have",
        xlabel="Acceleration factor", ylabel="Iterations to converge",
        logy=True,
        note="Missing points are non-convergence: IEEE 39-bus fails above 1.4.")
    print("  ", V.save(fig, "wp1_gs_acceleration"))


def fig_wp2_system_strength() -> None:
    rows = load("wp2_system_strength.csv")
    if not rows:
        return
    by_case = defaultdict(list)
    for r in rows:
        by_case[r["case"]].append(r)

    for case, sub in by_case.items():
        sub = sorted(sub, key=lambda r: fnum(r, "actual_pen_pct"))
        pens = [fnum(r, "actual_pen_pct") for r in sub]
        series = {
            "Minimum SCR": [fnum(r, "min_scr") for r in sub],
            "WSCR (ERCOT weighted)": [fnum(r, "wscr") for r in sub],
            "CSCR (composite)": [fnum(r, "cscr") for r in sub],
        }
        # drop the infinite 0% point, which has no IBR to rate
        keep = [i for i, p in enumerate(pens) if p > 0]
        pens = [pens[i] for i in keep]
        series = {k: [v[i] for i in keep] for k, v in series.items()}
        if not pens:
            continue

        fig = V.line_by_series(
            pens, series,
            title=f"System strength collapses with IBR penetration "
                  f"({case.replace('case', 'IEEE ')}-bus)",
            xlabel="IBR penetration (% of installed capacity)",
            ylabel="Short-circuit ratio",
            threshold=(3.0, "SCR = 3.0 interconnection screen"),
            note="Below the screen, grid-following converters risk PLL instability "
                 "and sub-synchronous control interaction.")
        print("  ", V.save(fig, f"wp2_system_strength_{case}"))


def fig_wp2_penetration_loading() -> None:
    rows = load("wp2_penetration_loading.csv")
    if not rows:
        return
    for case in sorted({r["case"] for r in rows}):
        sub = [r for r in rows if r["case"] == case]
        loads = sorted({fnum(r, "load_scale") for r in sub})
        pens = sorted({fnum(r, "requested_pen_pct") for r in sub})

        grid = np.full((len(loads), len(pens)), np.nan)
        for i, lam in enumerate(loads):
            for j, p in enumerate(pens):
                r = next((x for x in sub
                          if abs(fnum(x, "load_scale") - lam) < 1e-9
                          and abs(fnum(x, "requested_pen_pct") - p) < 1e-9), None)
                if r and r["classical_converged"] == "True" and r["ibr_converged"] == "True":
                    grid[i, j] = fnum(r, "max_dvm_pu")

        fig = V.heatmap(
            grid,
            [f"x{l:.2f}" for l in loads], [f"{int(p)}%" for p in pens],
            title=f"Classical vs IBR-aware power flow: error appears only under stress "
                  f"({case.replace('case', 'IEEE ')}-bus)",
            xlabel="IBR penetration", ylabel="Load scaling",
            cbar_label="max |dVm| (pu)",
            note="The 0% column is identically zero by construction -- that is the "
                 "experimental control. Grey cells: one or both models did not converge.")
        print("  ", V.save(fig, f"wp2_penetration_loading_{case}"))


def fig_wp3_fault_error() -> None:
    rows = load("wp3_fault_error.csv")
    if not rows:
        return
    for case in sorted({r["case"] for r in rows}):
        sub = [r for r in rows if r["case"] == case and r["converged"] == "True"]
        if not sub:
            continue
        pens = sorted({fnum(r, "pen_pct") for r in sub})
        series = {}
        for ft in ["3LG", "SLG", "LL", "LLG"]:
            vals = []
            for p in pens:
                v = [fnum(r, "error_pct") for r in sub
                     if r["fault_type"] == ft and abs(fnum(r, "pen_pct") - p) < 1e-6]
                v = [x for x in v if np.isfinite(x)]
                vals.append(float(np.mean(v)) if v else float("nan"))
            series[ft] = vals

        fig = V.line_by_series(
            pens, series,
            title=f"Fault current falls as inverters displace machines "
                  f"({case.replace('case', 'IEEE ')}-bus)",
            xlabel="IBR penetration (% of installed capacity)",
            ylabel="Fault current error vs classical (%)",
            note="Converged solves only. Negative = the classical study OVERSTATES "
                 "available fault current, which is the direction that matters for "
                 "protection settings.")
        print("  ", V.save(fig, f"wp3_fault_error_{case}"))


def fig_wp4_misoperation() -> None:
    rows = load("wp4_protection_misoperation.csv")
    if not rows:
        return
    for case in sorted({r["case"] for r in rows}):
        sub = [r for r in rows if r["case"] == case]
        pens = sorted({fnum(r, "pen_pct") for r in sub})
        series = {}
        for tag, label in [("no", "Pre-2800 inverter (K2 = 0)"),
                           ("yes", "IEEE 2800 neg-seq injection (K2 = 4)")]:
            vals = []
            for p in pens:
                s = [r for r in sub if r["ieee2800"] == tag
                     and abs(fnum(r, "pen_pct") - p) < 1e-6]
                if not s:
                    vals.append(float("nan"))
                    continue
                n_bad = sum(1 for r in s if int(fnum(r, "n_misoperations")) > 0)
                vals.append(n_bad / len(s) * 100.0)
            series[label] = vals

        fig = V.line_by_series(
            pens, series,
            title=f"Relays set from a classical study increasingly misoperate "
                  f"({case.replace('case', 'IEEE ')}-bus)",
            xlabel="IBR penetration (% of installed capacity)",
            ylabel="Relay cases misoperating (%)",
            note="Dominant modes: overcurrent element runs slow, then fails to pick "
                 "up at all. IEEE 2800 injection helps at moderate penetration but "
                 "not at high, where fault current magnitude is the binding limit.")
        print("  ", V.save(fig, f"wp4_misoperation_{case}"))


def fig_wp5_gfm_mitigation() -> None:
    rows = load("wp5_gfm_wellposedness.csv")
    if rows:
        for case in sorted({r["case"] for r in rows}):
            sub = [r for r in rows if r["case"] == case]
            pens = sorted({fnum(r, "pen_pct") for r in sub})
            shares = sorted({fnum(r, "gfm_share_pct") for r in sub})
            series = {}
            for p in pens:
                vals = []
                for s in shares:
                    r = next((x for x in sub if abs(fnum(x, "pen_pct") - p) < 1e-6
                              and abs(fnum(x, "gfm_share_pct") - s) < 1e-6), None)
                    vals.append(fnum(r, "converged_pct") if r else float("nan"))
                series[f"{p:.0f}% IBR penetration"] = vals

            fig = V.line_by_series(
                shares, series,
                title=f"Grid-forming control restores a solvable fault problem "
                      f"({case.replace('case', 'IEEE ')}-bus)",
                xlabel="Share of the IBR fleet that is grid-forming (% of MVA)",
                ylabel="Fault cases solved (%)",
                note="At 0% grid-forming and 100% IBR there is no voltage source "
                     "anywhere in the network, so the phasor fault problem is "
                     "ill-posed. 25% grid-forming is enough to restore it.")
            print("  ", V.save(fig, f"wp5_gfm_wellposedness_{case}"))

    rows = load("wp5_gfm_strength.csv")
    if rows:
        for case in sorted({r["case"] for r in rows}):
            sub = [r for r in rows if r["case"] == case]
            pens = sorted({fnum(r, "pen_pct") for r in sub})
            shares = sorted({fnum(r, "gfm_share_pct") for r in sub})
            series = {}
            for p in pens:
                vals = []
                for s in shares:
                    r = next((x for x in sub if abs(fnum(x, "pen_pct") - p) < 1e-6
                              and abs(fnum(x, "gfm_share_pct") - s) < 1e-6), None)
                    vals.append(fnum(r, "wscr") if r else float("nan"))
                series[f"{p:.0f}% IBR penetration"] = vals

            fig = V.line_by_series(
                shares, series,
                title=f"System strength recovers with grid-forming share "
                      f"({case.replace('case', 'IEEE ')}-bus)",
                xlabel="Share of the IBR fleet that is grid-forming (% of MVA)",
                ylabel="WSCR (ERCOT weighted short-circuit ratio)",
                threshold=(3.0, "SCR = 3.0 interconnection screen"),
                note="A grid-forming converter is a voltage source behind a virtual "
                     "impedance, so it contributes system strength the way a "
                     "synchronous machine does.")
            print("  ", V.save(fig, f"wp5_gfm_strength_{case}"))


def fig_wp5_syncon() -> None:
    rows = load("wp5_syncon_mitigation.csv")
    if not rows:
        return
    for case in sorted({r["case"] for r in rows}):
        sub = sorted([r for r in rows if r["case"] == case],
                     key=lambda r: fnum(r, "syncon_mva_each"))
        mva = [fnum(r, "syncon_mva_each") for r in sub]
        fig = V.line_by_series(
            mva,
            {"Short-circuit level at the weak sites (MVA)":
                [fnum(r, "ssc_at_sites_mva") for r in sub]},
            title=f"Synchronous condensers restore short-circuit strength "
                  f"({case.replace('case', 'IEEE ')}-bus, 80% IBR)",
            xlabel="Synchronous condenser rating installed at each of 2 sites (MVA)",
            ylabel="Short-circuit level (MVA)",
            note="Sited at the weakest buses that host no generator -- which is "
                 "where a utility would actually place them.")
        print("  ", V.save(fig, f"wp5_syncon_ssc_{case}"))

        fig = V.line_by_series(
            mva,
            {"RoCoF after N-1 loss of the largest unit":
                [abs(fnum(r, "rocof_hz_s")) for r in sub]},
            title=f"Condenser inertia cuts RoCoF "
                  f"({case.replace('case', 'IEEE ')}-bus, 80% IBR)",
            xlabel="Synchronous condenser rating installed at each of 2 sites (MVA)",
            ylabel="|RoCoF| (Hz/s)",
            threshold=(1.0, "1.0 Hz/s loss-of-mains protection limit"),
            note="A condenser has no prime mover, so it adds no energy -- but its "
                 "rotating mass adds inertia, which is what bounds RoCoF.")
        print("  ", V.save(fig, f"wp5_syncon_rocof_{case}"))


def fig_wp6_cpf() -> None:
    rows = load("wp6_cpf_margin.csv")
    if not rows:
        return
    cases = sorted({r["case"] for r in rows})
    pens = sorted({fnum(r, "requested_pen_pct") for r in rows})
    series = {}
    for c in cases:
        vals = []
        for p in pens:
            r = next((x for x in rows if x["case"] == c
                      and abs(fnum(x, "requested_pen_pct") - p) < 1e-6), None)
            vals.append(fnum(r, "loading_margin_pct") if r else float("nan"))
        series[c.replace("case", "IEEE ") + "-bus"] = vals

    fig = V.line_by_series(
        pens, series,
        title="Voltage stability margin vs IBR penetration",
        xlabel="IBR penetration (% of installed capacity)",
        ylabel="Loading margin to collapse (%)",
        note="Converter capability applied to displaced machines. At part load a "
             "converter's reactive capability exceeds a fixed nameplate Qmax, so "
             "the classical fixed-limit model is conservative here.")
    print("  ", V.save(fig, "wp6_cpf_margin"))


def fig_wp6_n1() -> None:
    """N-1 security: how many outages newly break a limit."""
    rows = load("wp6_n1_ranking.csv")
    if not rows:
        return
    cases = sorted({r["case"] for r in rows})
    labels, solved_clean, solved_viol, islanded, diverged = [], [], [], [], []

    for c in cases:
        for pen in sorted({fnum(r, "pen_pct") for r in rows if r["case"] == c}):
            sub = [r for r in rows if r["case"] == c and abs(fnum(r, "pen_pct") - pen) < 1e-6]
            solved = [r for r in sub if r["status"] == "solved"]
            viol = sum(1 for r in solved if int(fnum(r, "n_violations")) > 0)
            labels.append(f"{c.replace('case', '')}\n{pen:.0f}% IBR")
            solved_clean.append(len(solved) - viol)
            solved_viol.append(viol)
            islanded.append(sum(1 for r in sub if r["status"] == "islands network"))
            diverged.append(sum(1 for r in sub if r["status"] == "diverged"))

    fig = V.grouped_bars(
        labels,
        {"Solved, clean": solved_clean,
         "New violations": solved_viol,
         "Islanded": islanded,
         "Diverged": diverged},
        title="N-1 branch-outage outcomes",
        xlabel="Test system and IBR penetration",
        ylabel="Number of branch outages",
        note="Violations are counted relative to the intact case -- only limits the "
             "outage NEWLY breaks. Islanding is separated from divergence: it is a "
             "known topological fact, not a numerical failure.")
    print("  ", V.save(fig, "wp6_n1_outcomes"))


def fig_wp6_l_index() -> None:
    rows = load("wp6_l_index.csv")
    if not rows:
        return
    for case in sorted({r["case"] for r in rows}):
        sub = [r for r in rows if r["case"] == case]
        pens = sorted({fnum(r, "pen_pct") for r in sub})
        series = {}
        for lam in sorted({fnum(r, "load_scale") for r in sub}):
            vals = []
            for p in pens:
                r = next((x for x in sub if abs(fnum(x, "load_scale") - lam) < 1e-9
                          and abs(fnum(x, "pen_pct") - p) < 1e-6), None)
                vals.append(fnum(r, "max_l_index") if r else float("nan"))
            series[f"load x{lam:.2f}"] = vals
        if not pens:
            continue
        fig = V.line_by_series(
            pens, series,
            title=f"Worst-bus L-index vs IBR penetration "
                  f"({case.replace('case', 'IEEE ')}-bus)",
            xlabel="IBR penetration (% of installed capacity)",
            ylabel="L-index (0 = unloaded, 1 = collapse)",
            note="Computed from a single solved power flow, which is what makes the "
                 "L-index usable for online monitoring where a CPF is too slow.")
        print("  ", V.save(fig, f"wp6_l_index_{case}"))


def fig_wp7_virtual_inertia() -> None:
    """The controlled GFM result: CCT vs synthesised inertia."""
    rows = load("wp7_virtual_inertia_sweep.csv")
    if not rows:
        return
    series = {}
    h_vals = sorted({fnum(r, "h_virtual_s") for r in rows})
    for case in sorted({r["case"] for r in rows}):
        for pen in sorted({fnum(r, "pen_pct") for r in rows if r["case"] == case}):
            vals = []
            for h in h_vals:
                r = next((x for x in rows if x["case"] == case
                          and abs(fnum(x, "pen_pct") - pen) < 1e-6
                          and abs(fnum(x, "h_virtual_s") - h) < 1e-9), None)
                vals.append(fnum(r, "cct_s") if r else float("nan"))
            label = f"{case.replace('case', 'IEEE ')}-bus @ {pen:.0f}% IBR"
            series[label] = vals

    fig = V.line_by_series(
        h_vals, series,
        title="Grid-forming virtual inertia extends critical clearing time",
        xlabel="Grid-forming virtual inertia constant H (s)",
        ylabel="Critical clearing time (s)",
        note="CONTROLLED comparison: topology, dispatch and dynamic-unit count all "
             "held fixed, so only synthesised inertia varies. A grid-following "
             "penetration sweep is NOT a valid comparison here -- see wp7_cct_artefact.")
    print("  ", V.save(fig, "wp7_virtual_inertia_cct"))


def fig_wp7_artefact() -> None:
    """The negative result: why a naive GFL sweep misleads."""
    rows = load("wp7_cct_metric_artefact.csv")
    if not rows:
        return
    rows = sorted(rows, key=lambda r: fnum(r, "pen_pct"))
    pens = [fnum(r, "pen_pct") for r in rows]
    fig = V.grouped_bars(
        [f"{p:.0f}%" for p in pens],
        {"CCT (s)": [fnum(r, "cct_s") for r in rows],
         "Rotors remaining": [fnum(r, "n_dynamic_units") / 10.0 for r in rows]},
        title="Why a grid-following penetration sweep must NOT be read as a CCT result",
        xlabel="Grid-following IBR penetration",
        ylabel="CCT (s)  /  rotors remaining (tens)",
        value_fmt="{:.2f}",
        note="CCT appears to improve only because converting a machine to a "
             "grid-following converter deletes its rotor, and rotor-angle "
             "separation is the criterion being measured. Only the 0% bar is a "
             "valid CCT.")
    print("  ", V.save(fig, "wp7_cct_artefact"))


def fig_wp8_fdia() -> None:
    """Detectability of three attack classes vs magnitude."""
    rows = load("wp8_attack_detection.csv")
    if not rows:
        return
    for case in sorted({r["case"] for r in rows}):
        sub = [r for r in rows if r["case"] == case]
        thr = fnum(sub[0], "chi2_threshold")
        mags = sorted({fnum(r, "magnitude_deg") for r in sub
                       if r["attack"].startswith("fdia") and np.isfinite(fnum(r, "magnitude_deg"))})
        if not mags:
            continue
        series = {}
        for atk, label in [("fdia_linearised", "FDIA, linearised (DC model)"),
                           ("fdia_exact", "FDIA, exact (AC model)")]:
            vals = []
            for m in mags:
                r = next((x for x in sub if x["attack"] == atk
                          and abs(fnum(x, "magnitude_deg") - m) < 1e-9), None)
                vals.append(fnum(r, "objective_J") if r else float("nan"))
            series[label] = vals

        fig = V.line_by_series(
            mags, series,
            title=f"An exact AC false-data injection leaves no residual to detect "
                  f"({case.replace('case', 'IEEE ')}-bus)",
            xlabel="Intended state corruption (degrees of voltage angle)",
            ylabel="Weighted residual J(x)",
            threshold=(thr, "chi-squared detection threshold"),
            logy=True,
            note="The exact attack tracks the clean value at every magnitude, so the "
                 "residual test cannot fire at any threshold. The linearised (DC) "
                 "construction leaves a second-order residual and IS caught.")
        print("  ", V.save(fig, f"wp8_fdia_detectability_{case}"))


def main() -> None:
    config.ensure_output_dirs()
    print("Generating figures into", config.FIGURE_DIR)
    for fn in (fig_wp1_convergence, fig_wp1_gs_acceleration,
               fig_wp2_system_strength, fig_wp2_penetration_loading,
               fig_wp3_fault_error, fig_wp4_misoperation,
               fig_wp5_gfm_mitigation, fig_wp5_syncon,
               fig_wp6_cpf, fig_wp6_n1, fig_wp6_l_index,
               fig_wp7_virtual_inertia, fig_wp7_artefact, fig_wp8_fdia):
        try:
            fn()
        except Exception as exc:                       # noqa: BLE001
            print(f"  !! {fn.__name__}: {type(exc).__name__}: {exc}")
    print("done")


if __name__ == "__main__":
    main()
