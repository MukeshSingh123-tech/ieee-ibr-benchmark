"""Export scenarios.yaml / tolerances.yaml to JSON for the MATLAB side.

MATLAB's YAML support varies by release, but `jsondecode` has been stable
forever. Rather than maintain a second copy of the configuration in .m form --
which would drift, and would destroy the "single source of truth" property the
whole project rests on -- the YAML is mechanically translated to JSON here.

Run this whenever config/*.yaml changes, BEFORE running anything in MATLAB.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gridbench import config, metrics  # noqa: E402


def main() -> None:
    config.INTERCHANGE_DIR.mkdir(parents=True, exist_ok=True)

    targets = [
        ("scenarios.json", config.scenarios()),
        ("tolerances.json", config.tolerances()),
    ]
    for name, payload in targets:
        path = metrics.write_json(config.INTERCHANGE_DIR / name, payload)
        print(f"wrote {path}")

    print("\nMATLAB can now read these with:")
    print("    cfg = gb_config();          % scenarios")
    print("    tol = gb_config('tolerances');")


if __name__ == "__main__":
    main()
