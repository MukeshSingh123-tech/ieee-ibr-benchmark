"""Run every Python study and regenerate every table and figure, from scratch.

    cd IEEE/python
    python run_all.py              # everything
    python run_all.py --fast       # skip the slow studies (WP6, WP7)
    python run_all.py --list       # show what would run

The rule this enforces: if a number or a figure in the report cannot be produced
by this command, it does not go in the report.

Studies are ordered by dependency -- export_config first (it also generates the
MATLAB bridge), figures last (they read the committed CSVs rather than
recomputing, so the tables and the plots cannot drift apart).
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
import traceback
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")

# (module, description, slow?)
STUDIES = [
    ("studies.export_config", "export config to JSON for MATLAB", False),
    ("studies.wp1_baseline", "WP1  classical load flow + solver comparison", False),
    ("studies.wp2_penetration", "WP2  IBR power flow + system strength", False),
    ("studies.wp3_faults", "WP3/4 faults + protection misoperation", True),
    ("studies.wp5_mitigation", "WP5  GFM and synchronous condenser mitigation", True),
    ("studies.wp6_stability", "WP6  CPF, L-index, PTDF/LODF, N-1", True),
    ("studies.wp7_transient", "WP7  transient stability and CCT", True),
    ("studies.wp8_fdia", "WP8  state estimation and FDIA", True),
    ("studies.make_figures", "     regenerate all figures from the tables", False),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fast", action="store_true", help="skip the slow studies")
    ap.add_argument("--list", action="store_true", help="list studies and exit")
    ap.add_argument("--only", metavar="NAME",
                    help="run one study, e.g. --only wp3_faults")
    args = ap.parse_args()

    if args.list:
        for mod, desc, slow in STUDIES:
            print(f"  {'[slow] ' if slow else '       '}{mod.split('.')[-1]:<18} {desc}")
        return 0

    queue = STUDIES
    if args.only:
        queue = [s for s in STUDIES if s[0].endswith(args.only)]
        if not queue:
            print(f"no study matching {args.only!r}; try --list")
            return 2
    elif args.fast:
        queue = [s for s in STUDIES if not s[2]]

    print("=" * 78)
    print(f"Running {len(queue)} studies")
    print("=" * 78)

    failures = []
    t_start = time.perf_counter()

    for mod_name, desc, _slow in queue:
        print(f"\n{'-' * 78}\n>>> {desc}\n{'-' * 78}")
        t0 = time.perf_counter()
        try:
            module = importlib.import_module(mod_name)
            module.main()
            print(f"    [{time.perf_counter() - t0:.1f}s]")
        except Exception:                                   # noqa: BLE001
            failures.append(mod_name)
            print(f"    FAILED after {time.perf_counter() - t0:.1f}s")
            traceback.print_exc(limit=3)

    total = time.perf_counter() - t_start
    print("\n" + "=" * 78)
    print(f"{len(queue) - len(failures)}/{len(queue)} studies completed "
          f"in {total / 60:.1f} min")
    if failures:
        print("FAILED: " + ", ".join(failures))
    else:
        from gridbench import config
        n_tab = len(list(config.TABLE_DIR.glob("*.csv")))
        n_fig = len(list(config.FIGURE_DIR.glob("*.png")))
        print(f"{n_tab} tables in {config.TABLE_DIR}")
        print(f"{n_fig} figures in {config.FIGURE_DIR}")
        print("\nNext: python -m pytest tests/ -q")
    print("=" * 78)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
