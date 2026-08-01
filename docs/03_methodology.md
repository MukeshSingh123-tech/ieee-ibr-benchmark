# Methodology and assumptions

Every model in this project, what it assumes, and where it stops being valid.
Referenced from the code (`gridbench/inertia.py`, `gridbench/transient.py`) so a
reader who hits a limitation in the source can find its justification here.

**The organising principle:** the IEEE test cases carry only positive-sequence
load-flow data. Everything beyond that — machine reactances, zero-sequence
impedances, inertia constants, converter ratings — is either *published benchmark
data* (stated, cited) or an *assumption* (stated, declared in
`config/scenarios.yaml`, never buried in code). Which one it is, is recorded
below.

---

## 1. Data provenance

| Quantity | Source | Status |
|---|---|---|
| Bus/branch/gen load-flow data | `pandapower.networks` (MATPOWER cases) | published |
| Machine H, X′d for case9 | Anderson & Fouad, WSCC 3-machine | published |
| Machine H, X′d for case39 | Athay et al., New England 10-machine | published |
| Subtransient X″d, X₂, X₀ | typical utility values (IEEE C37.010) | **assumed** |
| Zero-sequence line multipliers | R₀/R₁ = 3.0, X₀/X₁ = 3.0, B₀/B₁ = 0.5 | **assumed** |
| Transformer winding connections | all delta-wye(g) → block zero sequence | **assumed** |
| Machine MVA rating | inferred from PMAX (pandapower leaves MBASE = NaN) | **inferred** |
| Slack short-circuit level | 10 × total system load | **assumed** |
| Converter current limit | 1.20 pu transient, 1.00 pu continuous | **assumed** |
| IEEE 2800 K₂ | 4.0, swept 0–6 | **assumed** |

The assumed values are all swept somewhere in the studies, so their influence on
the conclusions is visible rather than hidden.

---

## 2. Power flow

Standard polar Newton-Raphson, Gauss-Seidel, and fast-decoupled (XB/BX), written
from first principles in `gridbench/solvers.py` and `matlab/matpower/`.

**Assumptions:** balanced three-phase, positive sequence only, constant-power
loads, no frequency dependence.

**Validation:** hand-written solvers match pandapower to < 1e-8 pu; Ybus matches
to 5.7e-14. Gauss-Seidel is gated at 1e-6 pu because it converges on |ΔV| to
1e-8 rather than on power mismatch — holding it to the Newton gate would fail it
for being asked a different question.

---

## 3. IBR-aware power flow (WP2)

The **only** modelling change from classical is where the reactive limit sits:

```
classical:   Qmax = constant, from the case file
IBR-aware:   Qmax(V, P) = sqrt((V · Ilim · S)² − P²)
```

**The experimental control.** The replacement inverter is sized
`S = sqrt(PMAX² + QMAX²)` so that at V = 1.0 and P = PMAX its capability equals
the displaced machine's case-file QMAX *exactly* (verified to 1e-14 by
`test_inverter_sizing_matches_machine_capability`). Without this, an
independently inferred nameplate gives the inverter more capability than the
machine it replaces — it did, 33 vs 24 MVAr on IEEE 14-bus buses 6 and 8 — and
the study silently measures the sizing assumption instead of the physics.

**Consequence:** at 0% penetration the two models are the same computation, bit
for bit. Any nonzero difference is attributable to the limit's voltage and power
dependence and to nothing else.

**Limitation:** steady-state only. Converter dynamics, PLL behaviour and
switching are out of scope by construction.

---

## 4. Fault analysis (WP3)

Symmetrical components with the sequence networks built in
`gridbench/faults.py`. Three source types, treated as what they physically are:

| Source | Fault-network representation |
|---|---|
| Synchronous machine | voltage source behind X″d → shunt + Norton source |
| Grid-forming converter | voltage source behind Z_virtual → shunt + Norton source |
| Grid-following converter | current source, **no** shunt, no internal EMF |

Both solvers share one boundary-condition routine written in terms of
open-circuit sequence voltages at the fault bus. The classical case is the
special case V₂,oc = V₀,oc = 0 — which is what lets the two be compared without
a second implementation.

**The nonlinearity.** A grid-following converter's injection depends on the very
terminal voltage the fault produces, and is hard-clamped at the current limit. So
the network is nonlinear and must be iterated; superposition does not hold.

**PLL coast-through.** Below ~0.1 pu terminal voltage the injection angle
reference V/|V| is numerically meaningless (voltages reach 4e-5 pu at a close-in
bolted fault) and jitters, so no fixed point exists. Real inverters lose PLL lock
and coast on the last angle; modelling that is both physically correct and what
makes the problem well-conditioned. Blended over 0.02–0.10 pu to stay continuous.

**Well-posedness.** With no voltage source anywhere — no synchronous plant and no
grid-forming converters — every source is a current injection whose magnitude
depends on a voltage nothing establishes. This is reported as ill-posed rather
than iterated to a meaningless number.

**Limitations:** fundamental frequency only; no DC offset, no harmonics, no
sub-cycle behaviour, no converter trip. Zero-sequence data is synthesised.

---

## 5. Protection (WP4)

Textbook elements — mho distance (3 zones), negative-sequence directional (32Q),
IEC inverse-time overcurrent. Deliberately standard: the point is not a novel
relay, it is that **correctly-set conventional relays misoperate** when the
source behind them stops being a machine.

Relays are commissioned from the *classical* study, then evaluated against the
*IBR-aware* fault solution — which is exactly what happens in practice.

**Limitation:** steady-state phasor evaluation of each element. No CT saturation,
no relay filtering or algorithm latency, no pilot schemes. The dynamic sequence
(pickup → time → trip → reclose → lockout) lives in the Stateflow chart instead.

---

## 6. Frequency response (WP5)

Low-order System Frequency Response model: single centre-of-inertia frequency,
aggregated governor with reheat, load damping.

RoCoF is exact from the swing equation. The nadir comes from **numerical
integration** of the SFR ODE, not the Anderson-Mirheydar closed form — that
closed form assumes a moderately damped second-order response and returns
nonsense as H → 0, reporting the nadir *improving* at 84% penetration, which is
backwards. Numerical integration stays valid across the whole sweep.

**Limitations:** no network dynamics, no voltage coupling, no unit-by-unit
governor differences, no fast frequency response from converters.

---

## 7. Transient stability (WP7)

Classical multi-machine model: constant E′ behind X′d, swing equation per unit,
Kron-reduced to the internal source nodes. Grid-forming converters obey the same
swing equation with a synthesised inertia constant; grid-following converters
have no rotor and fold into the network as constant-impedance injections.

**The critical limitation, stated in the code as well as here:**

> **CCT is not comparable across scenarios with different numbers of dynamic
> units.** The criterion is rotor-angle separation, defined only over units that
> *have* a rotor angle. Converting a machine to a grid-following converter
> deletes it from the swing model, so the metric loses the failure mode it exists
> to detect. Measured on IEEE 39-bus: CCT "improves" from 0.121 s at 0% IBR to
> 0.977 s at 40.6% grid-following penetration. That is an artefact.

A grid-following converter destabilises through **PLL loss of synchronisation** —
a converter-control phenomenon on a millisecond timescale that an
electromechanical model cannot represent at all. This is an independent route to
the same conclusion the fault work reached: EMT analysis is mandatory in
converter-dominated grids.

**The valid comparison** is `virtual_inertia_sweep`: fixed topology, fixed
dispatch, fixed unit count, only the grid-forming inertia varies.

**Other limitations:** no exciter, no governor, no saliency, no damper windings,
loads as constant impedance. A screening model, validated against a published
CCT (0.121 s vs ~0.15 s on IEEE 39-bus) but not a substitute for EMT.

---

## 8. State estimation and FDIA (WP8)

Weighted least squares by Gauss-Newton, numerical measurement Jacobian,
chi-squared bad-data detection at α = 0.01, normalised residuals for
identification.

**Attack constructions:**

- *Random* — arbitrary perturbation of a few measurements. The baseline a
  defender should catch, and does.
- *Linearised*, `a = Hc` — the classic DC-model construction. **Not** stealthy
  against an AC estimator: `h(x+c) ≠ h(x) + Hc`, and the second-order remainder
  leaves a residual that grows with |c|. Caught from 2° upward.
- *Exact*, `a = h(x+c) − h(x)` — every measurement moved to precisely what it
  would read if the state really were `x + c`. Residual **unchanged**;
  undetectable at any magnitude or threshold.

**The claim being made** is narrow and structural: residual-based detection
cannot see an attack that lies in the range of the measurement function, because
there is no residual to test. It is not a claim that such attacks are easy to
mount — they require accurate topology and parameter knowledge, and the
linearised result shows what happens to an attacker who lacks it.

**Limitations:** static estimation (no forecasting or temporal correlation), no
PMU measurements, no topology-error processing, and no attempt at an actual
detector — the point is to quantify the gap that out-of-band verification exists
to close.

---

## 9. What is deliberately out of scope

| Not modelled | Why, and where it would go |
|---|---|
| EMT / switching transients | needs PSCAD or Simscape at 50 µs — deferred, scaffolded in `config` |
| PLL small-signal stability | the actual GFL failure mode; requires impedance-based analysis |
| Sub-synchronous control interaction | EMT or frequency-domain impedance scanning |
| Unbalanced distribution | OpenDSS, IEEE 13/123-node — deferred |
| Converter trip on ride-through failure | Stateflow FRT chart has the logic; not coupled to the phasor studies |
| Economic dispatch / LMP | OPF layer, not yet built |
| Harmonics | fundamental frequency only throughout |

---

## 10. Reproducibility

- Every scenario parameter lives in `config/scenarios.yaml`; every threshold in
  `config/tolerances.yaml`. Nothing is hard-coded in a study script.
- MATLAB reads the same YAML through a generated JSON bridge, so the two
  toolchains cannot drift apart.
- Monte Carlo seeds are fixed (`monte_carlo.seed: 20260801`).
- `python run_all.py` regenerates every table and figure from scratch.
- Figures are built from the committed CSVs, not recomputed — so a number in the
  report and the same number in a plot cannot disagree.
- Hypotheses were recorded in `tolerances.yaml` *before* running, with outcomes
  (CONFIRMED / REFUTED) written back afterwards.
