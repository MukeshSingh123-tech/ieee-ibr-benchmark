"""Bus admittance and impedance matrix construction, written from first principles.

This is deliberately NOT delegated to pandapower or MATPOWER. Ybus and Zbus are
the foundation of every study downstream (load flow, fault, system strength),
and `validate_against_reference()` proves this implementation reproduces
pandapower's internal Ybus to machine precision -- so the "from scratch" claim
is checkable, not asserted.

Branch model (MATPOWER convention), for a line from bus f to bus t with series
impedance r + jx, total line charging susceptance b, off-nominal tap ratio tau
and phase shift theta (degrees), with the tap on the FROM side:

    ys = 1 / (r + jx)          series admittance
    t  = tau * exp(j*theta)    complex tap

    Yff = (ys + j*b/2) / |t|^2
    Yft = -ys / conj(t)
    Ytf = -ys / t
    Ytt =  ys + j*b/2
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .ppc import (
    BR_B, BR_R, BR_STATUS, BR_X, BS, BUS_I, F_BUS, GS, SHIFT, T_BUS, TAP,
    bus_index,
)


def branch_admittances(ppc: dict[str, Any]) -> tuple[np.ndarray, ...]:
    """Per-branch (Yff, Yft, Ytf, Ytt) primitive admittances."""
    br = ppc["branch"]
    status = br[:, BR_STATUS] > 0

    ys = np.zeros(br.shape[0], dtype=complex)
    z = br[:, BR_R] + 1j * br[:, BR_X]
    # out-of-service branches contribute nothing; guard against z == 0
    nz = status & (np.abs(z) > 0)
    ys[nz] = 1.0 / z[nz]

    bc = np.where(status, br[:, BR_B], 0.0)

    tau = br[:, TAP].copy()
    tau[tau == 0] = 1.0            # MATPOWER: TAP == 0 means nominal ratio
    shift = np.deg2rad(br[:, SHIFT])
    t = tau * np.exp(1j * shift)

    yff = (ys + 1j * bc / 2.0) / (tau ** 2)
    yft = -ys / np.conj(t)
    ytf = -ys / t
    ytt = ys + 1j * bc / 2.0
    return yff, yft, ytf, ytt


def build(ppc: dict[str, Any]) -> np.ndarray:
    """Dense complex bus admittance matrix in per-unit on system base."""
    bus, br = ppc["bus"], ppc["branch"]
    n = bus.shape[0]
    idx = bus_index(ppc)

    ybus = np.zeros((n, n), dtype=complex)
    yff, yft, ytf, ytt = branch_admittances(ppc)

    for k in range(br.shape[0]):
        f = idx[int(br[k, F_BUS])]
        t = idx[int(br[k, T_BUS])]
        ybus[f, f] += yff[k]
        ybus[f, t] += yft[k]
        ybus[t, f] += ytf[k]
        ybus[t, t] += ytt[k]

    # fixed bus shunts, given in MW / MVAr demanded at V = 1.0 pu
    ybus[np.arange(n), np.arange(n)] += (bus[:, GS] + 1j * bus[:, BS]) / ppc["baseMVA"]
    return ybus


def build_zbus(ppc: dict[str, Any], ybus: np.ndarray | None = None) -> np.ndarray:
    """Zbus = inv(Ybus).

    Direct inversion is fine at IEEE-test-system scale (<= 123 buses) and is far
    clearer than the sequential Zbus building algorithm. For fault studies Zbus
    must be a full matrix anyway, so sparsity buys nothing here.
    """
    ybus = build(ppc) if ybus is None else ybus
    n = ybus.shape[0]
    if np.linalg.matrix_rank(ybus) < n:
        raise np.linalg.LinAlgError(
            "Ybus is singular -- the network has no shunt path to ground. "
            "Add a shunt/ground reference before forming Zbus."
        )
    return np.linalg.inv(ybus)


def validate_against_reference(ppc: dict[str, Any], atol: float = 1e-9) -> float:
    """Compare our Ybus with pandapower's internal one. Returns max abs error.

    Raises AssertionError if the mismatch exceeds `atol`. This is the gate that
    makes every downstream result trustworthy: if Ybus is wrong, everything is.
    """
    ref = ppc.get("_ref_Ybus")
    if ref is None:
        raise KeyError("ppc has no _ref_Ybus; build it with ppc.load_ppc()")
    err = float(np.max(np.abs(build(ppc) - ref)))
    assert err <= atol, f"Ybus mismatch {err:.3e} exceeds {atol:.1e}"
    return err


def thevenin_impedance(ppc: dict[str, Any], bus: int, zbus: np.ndarray | None = None) -> complex:
    """Thevenin (driving-point) impedance looking into `bus`, in per-unit.

    This is Zbus[i, i]. It is the quantity PSCAD needs for the reduced EMT
    equivalent (WP5) and the quantity SCR is computed from (WP2), so the phasor
    and EMT halves of the project are anchored to the same number.
    """
    zbus = build_zbus(ppc) if zbus is None else zbus
    i = bus_index(ppc)[int(bus)]
    return complex(zbus[i, i])


def branch_flows(ppc: dict[str, Any], v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Complex power flow (Sf, St) at both ends of every branch, per-unit.

    Sign convention: positive = flowing INTO the branch at that end, so
    Sf + St equals the branch losses.
    """
    br = ppc["branch"]
    idx = bus_index(ppc)
    f = np.array([idx[int(x)] for x in br[:, F_BUS]])
    t = np.array([idx[int(x)] for x in br[:, T_BUS]])

    yff, yft, ytf, ytt = branch_admittances(ppc)
    i_f = yff * v[f] + yft * v[t]
    i_t = ytf * v[f] + ytt * v[t]
    return v[f] * np.conj(i_f), v[t] * np.conj(i_t)


def total_losses_mw(ppc: dict[str, Any], v: np.ndarray) -> float:
    sf, st = branch_flows(ppc, v)
    return float(np.real(sf + st).sum() * ppc["baseMVA"])
