"""Error metrics, angle handling, and table output shared by every study.

Kept in one place so that "percentage error" and "angle difference" mean the
same thing in every table and figure in the report -- including the ones the
MATLAB side produces, which write the same column names.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


# =============================================================================
# angles
# =============================================================================
def wrap_deg(angle: float | np.ndarray) -> np.ndarray | float:
    """Wrap an angle (or array) into (-180, +180] degrees.

    Necessary for every angle comparison in this project: a directional relay
    that sees -161.4 deg and 177.0 deg differs by -21.5 deg, not +338.5 deg.
    Reporting the unwrapped value would turn a modest, plausible shift into a
    spectacular but meaningless one.
    """
    wrapped = (np.asarray(angle, dtype=float) + 180.0) % 360.0 - 180.0
    # (-180, 180]: map the -180 endpoint to +180
    wrapped = np.where(np.isclose(wrapped, -180.0), 180.0, wrapped)
    return float(wrapped) if np.ndim(angle) == 0 else wrapped


def angle_diff_deg(a: float | np.ndarray, b: float | np.ndarray):
    """Shortest signed angular difference a - b, in degrees."""
    return wrap_deg(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))


# =============================================================================
# error metrics
# =============================================================================
def pct_error(actual: float | np.ndarray, reference: float | np.ndarray):
    """Signed percentage error of `actual` relative to `reference`.

    Returns NaN where the reference is ~0 rather than a huge number, so a
    near-zero denominator cannot masquerade as a dramatic finding.
    """
    ref = np.asarray(reference, dtype=float)
    act = np.asarray(actual, dtype=float)
    out = np.where(np.abs(ref) > 1e-12, (act - ref) / np.where(np.abs(ref) > 1e-12, ref, 1.0) * 100.0, np.nan)
    return float(out) if np.ndim(reference) == 0 and np.ndim(actual) == 0 else out


def max_abs_error(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b)))) if np.size(a) else float("nan")


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(d ** 2))) if d.size else float("nan")


def compare_voltages(v_a: np.ndarray, v_b: np.ndarray) -> dict[str, float]:
    """Standard voltage comparison block used by every cross-model table."""
    vm_a, vm_b = np.abs(v_a), np.abs(v_b)
    va_a, va_b = np.degrees(np.angle(v_a)), np.degrees(np.angle(v_b))
    return {
        "max_dvm_pu": max_abs_error(vm_a, vm_b),
        "rmse_vm_pu": rmse(vm_a, vm_b),
        "max_dva_deg": float(np.max(np.abs(angle_diff_deg(va_a, va_b)))),
        "vm_min_pu": float(np.min(vm_b)),
        "vm_max_pu": float(np.max(vm_b)),
    }


# =============================================================================
# output
# =============================================================================
def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], float_fmt: str = "%.6g") -> Path:
    """Write a list of dicts to CSV, creating parent directories as needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return path

    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                k: (float_fmt % v if isinstance(v, float) and np.isfinite(v) else v)
                for k, v in row.items()
            })
    return path


def write_json(path: Path, payload: Any) -> Path:
    """Write JSON, converting numpy scalars/arrays so the MATLAB side can read it."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def _default(obj: Any) -> Any:
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return None if not np.isfinite(obj) else float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, complex):
            return {"re": obj.real, "im": obj.imag}
        raise TypeError(f"cannot serialise {type(obj).__name__}")

    path.write_text(json.dumps(payload, indent=2, default=_default), encoding="utf-8")
    return path


def markdown_table(rows: Sequence[Mapping[str, Any]], fmt: str = "{:.4g}") -> str:
    """Render rows as a GitHub-flavoured markdown table for the README/report."""
    if not rows:
        return "_(no rows)_"
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    def cell(value: Any) -> str:
        if isinstance(value, float):
            return "n/a" if not np.isfinite(value) else fmt.format(value)
        return str(value)

    lines = ["| " + " | ".join(fields) + " |",
             "|" + "|".join("---" for _ in fields) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(cell(row.get(f, "")) for f in fields) + " |")
    return "\n".join(lines)


def print_table(rows: Sequence[Mapping[str, Any]], fmt: str = "{:>12.5g}") -> None:
    """Console table for interactive runs."""
    if not rows:
        print("(no rows)")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    print("  ".join(f"{f:>12}" for f in fields))
    for row in rows:
        cells = []
        for f in fields:
            v = row.get(f, "")
            if isinstance(v, float):
                cells.append("         n/a" if not np.isfinite(v) else fmt.format(v))
            else:
                cells.append(f"{str(v):>12}")
        print("  ".join(cells))
