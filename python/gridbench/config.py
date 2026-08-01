"""Project configuration: paths, scenarios.yaml, tolerances.yaml.

Every study script starts here. Nothing in this project may hard-code a
penetration level, fault location, current limit or tolerance -- it comes from
config/scenarios.yaml and config/tolerances.yaml or it does not exist.
"""

from __future__ import annotations

import functools
from pathlib import Path
from typing import Any

import yaml

# IEEE/python/gridbench/config.py -> IEEE/
ROOT = Path(__file__).resolve().parents[2]

CONFIG_DIR = ROOT / "config"
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROFILE_DIR = DATA_DIR / "profiles"
INTERCHANGE_DIR = DATA_DIR / "interchange"
RESULTS_DIR = ROOT / "results"
FIGURE_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"
LOG_DIR = RESULTS_DIR / "logs"
MATLAB_DIR = ROOT / "matlab"


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


@functools.lru_cache(maxsize=1)
def scenarios() -> dict[str, Any]:
    """Parsed config/scenarios.yaml (cached)."""
    return _read_yaml(CONFIG_DIR / "scenarios.yaml")


@functools.lru_cache(maxsize=1)
def tolerances() -> dict[str, Any]:
    """Parsed config/tolerances.yaml (cached)."""
    return _read_yaml(CONFIG_DIR / "tolerances.yaml")


def tol(*keys: str) -> float:
    """Fetch a nested tolerance, e.g. tol('cross_tool', 'matpower_vs_pandapower', 'vm_pu').

    Raises rather than defaulting: a missing tolerance is a config bug, and
    silently substituting one would let a study pass a gate that does not exist.
    """
    node: Any = tolerances()
    trail: list[str] = []
    for key in keys:
        trail.append(key)
        if not isinstance(node, dict) or key not in node:
            raise KeyError(f"tolerances.yaml has no entry {'.'.join(trail)}")
        node = node[key]
    if not isinstance(node, (int, float)):
        raise TypeError(f"tolerance {'.'.join(keys)} is {type(node).__name__}, not a number")
    return float(node)


def system(name: str) -> dict[str, Any]:
    """Definition block for one test system from scenarios.yaml."""
    systems = scenarios()["systems"]
    if name not in systems:
        raise KeyError(f"unknown system {name!r}; known: {sorted(systems)}")
    return systems[name]


def systems_for(tool: str) -> list[str]:
    """Which test systems a given tool can hold, honouring node-count caps."""
    return [
        name
        for name, spec in scenarios()["systems"].items()
        if tool in spec.get("tools", [])
    ]


def ensure_output_dirs() -> None:
    """Create results/ subtrees. Safe to call repeatedly."""
    for path in (FIGURE_DIR, TABLE_DIR, LOG_DIR, INTERCHANGE_DIR):
        path.mkdir(parents=True, exist_ok=True)


__all__ = [
    "ROOT", "CONFIG_DIR", "DATA_DIR", "RAW_DIR", "PROFILE_DIR", "INTERCHANGE_DIR",
    "RESULTS_DIR", "FIGURE_DIR", "TABLE_DIR", "LOG_DIR", "MATLAB_DIR",
    "scenarios", "tolerances", "tol", "system", "systems_for", "ensure_output_dirs",
]
