"""WP5 -- Mitigation: what actually fixes a converter-dominated grid.

WP1-WP4 quantify the problem. This study quantifies the two remedies a system
operator can actually buy, and shows that both work.

  1. GRID-FORMING CONTROL
     A grid-following converter is a current source: it needs a voltage to
     synchronise to. A grid-forming converter is a VOLTAGE source behind a
     virtual impedance -- it establishes the reference rather than tracking it.
     At 100% IBR with an all-grid-following fleet the phasor fault problem has
     no voltage source at all and is genuinely ill-posed. Converting part of the
     fleet to grid-forming restores it. That is not a numerical trick; it is the
     physical reason grid codes are beginning to mandate GFM capability.

  2. SYNCHRONOUS CONDENSERS
     A machine with no prime mover: no energy, but inertia, fault current and a
     voltage reference. This is what National Grid, EirGrid and AEMO have
     actually procured when system strength ran short.

Outputs
    results/tables/wp5_gfm_wellposedness.csv
    results/tables/wp5_gfm_strength.csv
    results/tables/wp5_syncon_mitigation.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from gridbench import (  # noqa: E402
    config, faults as F, ibr as I, inertia as IN, metrics, ppc as P,
    solvers as S, strength as ST,
)

CASES = ["case14", "case39"]
GFM_SHARES = [0, 25, 50, 75, 100]
FAULT_TYPES = ["3LG", "SLG", "LL", "LLG"]


def _networks(c, base, ibr_gens, gfm_share):
    gfl_g, gfm_g = I.split_gfl_gfm(c, ibr_gens, gfm_share)
    return F.build_sequence_networks(
        c, base.v,
        ibr_buses=I.ibr_buses_from_gens(c, ibr_gens),
        gfm_buses=I.ibr_buses_from_gens(c, gfm_g) if gfm_g.size else None,
    )


def gfm_wellposedness() -> list[dict]:
    """Does grid-forming control restore a solvable fault problem?"""
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        base = S.newton_raphson(c)
        locations = config.scenarios()["faults"]["locations"].get(case, [])
        for pen in [60, 80, 100]:
            gens = I.select_ibr_gens(c, pen)
            if gens.size == 0:
                continue
            for share in GFM_SHARES:
                seq = _networks(c, base, gens, share)
                n_ok = n_wp = n_tot = 0
                mags = []
                for bus in locations:
                    for ft in FAULT_TYPES:
                        r = F.solve_fault_ibr_aware(c, seq, bus, ft)
                        n_tot += 1
                        n_ok += bool(r.converged)
                        n_wp += bool(r.well_posed)
                        if r.converged:
                            mags.append(r.i_fault_mag)
                rows.append({
                    "case": case,
                    "pen_pct": round(I.actual_penetration_pct(c, gens), 1),
                    "gfm_share_pct": share,
                    "voltage_source_buses": int(seq.voltage_source_buses.size),
                    "n_gfm_buses": int(seq.gfm_buses.size),
                    "n_gfl_buses": int(seq.gfl_buses.size),
                    "well_posed_pct": n_wp / n_tot * 100 if n_tot else float("nan"),
                    "converged_pct": n_ok / n_tot * 100 if n_tot else float("nan"),
                    "mean_fault_current_pu": float(np.mean(mags)) if mags else float("nan"),
                })
    return rows


def gfm_strength() -> list[dict]:
    """Short-circuit level and SCR as grid-forming share increases."""
    rows = []
    for case in CASES:
        c = P.load_ppc(case)
        base = S.newton_raphson(c)
        for pen in [60, 80, 100]:
            gens = I.select_ibr_gens(c, pen)
            if gens.size == 0:
                continue
            for share in GFM_SHARES:
                seq = _networks(c, base, gens, share)
                prof = ST.profile(c, seq, gens)
                rows.append({
                    "case": case,
                    "pen_pct": round(I.actual_penetration_pct(c, gens), 1),
                    "gfm_share_pct": share,
                    "min_scr": prof["min_scr"],
                    "wscr": prof["wscr"],
                    "cscr": prof["cscr"],
                    "n_weak_buses": prof["n_weak_buses"],
                })
    return rows


def synchronous_condenser_mitigation() -> list[dict]:
    """Install synchronous condensers at weak NON-GENERATOR buses and re-measure.

    Two things this gets right, both of which were wrong on the first attempt and
    silently produced a null result:

      1. SITING. Condensers go at buses that do NOT already host an IBR. A bus
         hosting an IBR is classified as an IBR bus by the sequence-network
         builder, so a condenser placed there is excluded from the machine
         shunts and contributes nothing at all. Utilities site condensers at
         substations anyway, which is what this now models.

      2. A FIXED IBR FLEET. The IBR generators are chosen ONCE on the base case
         and reused. Re-selecting after adding condensers changes total system
         capacity and therefore changes which machines are displaced, so the
         comparison would no longer be like-for-like. Condensers are appended to
         the gen table, leaving the original row indices valid.
    """
    rows = []
    for case in CASES:
        c0 = P.load_ppc(case)
        base0 = S.newton_raphson(c0)
        gens0 = I.select_ibr_gens(c0, 80)          # chosen ONCE, reused throughout
        if gens0.size == 0:
            continue

        ibr_bus_numbers = {int(c0["gen"][g, P.GEN_BUS]) for g in gens0}
        gen_bus_numbers = {int(b) for b in c0["gen"][:, P.GEN_BUS]}
        seq0 = F.build_sequence_networks(
            c0, base0.v, ibr_buses=I.ibr_buses_from_gens(c0, gens0))

        # weakest buses that host neither a generator nor an IBR
        candidates = [
            int(b) for b in c0["bus"][:, P.BUS_I].astype(int)
            if b not in gen_bus_numbers and b not in ibr_bus_numbers
        ]
        candidates.sort(key=lambda b: F.short_circuit_mva(c0, seq0, b))
        targets = candidates[:2]
        if not targets:
            continue

        for mva in [0, 100, 200, 400]:
            c = c0 if mva == 0 else I.add_synchronous_condensers(c0, targets, mva)
            base = S.newton_raphson(c, enforce_q_limits=True)
            if not base.converged:
                continue
            seq = F.build_sequence_networks(
                c, base.v, ibr_buses=I.ibr_buses_from_gens(c, gens0))
            prof = ST.profile(c, seq, gens0)

            dp = IN.largest_single_contingency_mw(c)
            fr = IN.frequency_response(c, dp, gens0)

            # short-circuit level at the condenser sites, the direct measure
            ssc = [F.short_circuit_mva(c, seq, b) for b in targets]

            rows.append({
                "case": case,
                "syncon_mva_each": mva,
                "syncon_buses": " ".join(str(b) for b in targets),
                "ssc_at_sites_mva": float(np.mean(ssc)),
                "min_scr": prof["min_scr"],
                "wscr": prof["wscr"],
                "n_weak_buses": prof["n_weak_buses"],
                "h_sys_s": fr.h_sys_s,
                "rocof_hz_s": fr.rocof_hz_s,
                "violates_rocof": fr.violates_rocof,
            })
    return rows


def main() -> None:
    config.ensure_output_dirs()
    print("=" * 78)
    print("WP5 -- MITIGATION: grid-forming control and synchronous condensers")
    print("=" * 78)

    wp = gfm_wellposedness()
    metrics.write_csv(config.TABLE_DIR / "wp5_gfm_wellposedness.csv", wp)

    print("\nDoes grid-forming control restore a solvable fault problem?")
    print("  (share of the IBR fleet that is grid-forming, by MVA)")
    for case in CASES:
        pens = sorted({r["pen_pct"] for r in wp if r["case"] == case})
        for pen in pens:
            print(f"\n  {case} @ {pen}% IBR penetration")
            print(f"    {'GFM share':>10} {'V-source buses':>15} {'well-posed':>12} "
                  f"{'converged':>11} {'mean |If| pu':>13}")
            for share in GFM_SHARES:
                r = next((x for x in wp if x["case"] == case
                          and x["pen_pct"] == pen and x["gfm_share_pct"] == share), None)
                if r is None:
                    continue
                mag = (f"{r['mean_fault_current_pu']:13.4f}"
                       if np.isfinite(r["mean_fault_current_pu"]) else f"{'n/a':>13}")
                print(f"    {share:9d}% {r['voltage_source_buses']:15d} "
                      f"{r['well_posed_pct']:11.0f}% {r['converged_pct']:10.0f}%{mag}")

    st = gfm_strength()
    metrics.write_csv(config.TABLE_DIR / "wp5_gfm_strength.csv", st)
    print("\n\nSystem strength recovery with grid-forming share:")
    for case in CASES:
        pens = sorted({r["pen_pct"] for r in st if r["case"] == case})
        for pen in pens:
            cells = []
            for share in GFM_SHARES:
                r = next((x for x in st if x["case"] == case
                          and x["pen_pct"] == pen and x["gfm_share_pct"] == share), None)
                cells.append(f"{r['wscr']:8.2f}" if r else f"{'n/a':>8}")
            print(f"  {case:>8} @ {pen:5.1f}% IBR   WSCR: " +
                  "  ".join(f"{s}%:{c}" for s, c in zip(GFM_SHARES, cells)))

    sc = synchronous_condenser_mitigation()
    metrics.write_csv(config.TABLE_DIR / "wp5_syncon_mitigation.csv", sc)
    print("\n\nSynchronous condenser mitigation (at 80% IBR penetration):")
    print(f"  {'case':>8} {'MVA each':>9} {'Ssc@site':>10} {'min SCR':>9} {'WSCR':>8} "
          f"{'weak':>6} {'H_sys':>7} {'RoCoF':>8}")
    for r in sc:
        print(f"  {r['case']:>8} {r['syncon_mva_each']:9.0f} {r['ssc_at_sites_mva']:10.0f} "
              f"{r['min_scr']:9.2f} {r['wscr']:8.2f} {r['n_weak_buses']:6d} "
              f"{r['h_sys_s']:7.2f} {r['rocof_hz_s']:8.3f}")

    print(f"\nTables written to {config.TABLE_DIR}")


if __name__ == "__main__":
    main()
