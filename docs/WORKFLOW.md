# Workflow — what to actually do, in order

Two tracks. The **Python track** runs here in VS Code and produces the result
tables. The **MATLAB track** you run in the MATLAB app; it reproduces the same
studies independently and adds the dynamic work Python cannot do.

They meet at `results/tables/`, where every study writes a `*_matlab.csv`
alongside its Python `*.csv`, and the comparison harness diffs them.

---

## Track A — Python (VS Code)

### A1. One-time setup

```bash
cd IEEE/python
python -m pip install -r ../requirements.txt
python studies/export_config.py      # also generates the MATLAB config bridge
```

### A2. Run the studies

```bash
python studies/wp1_baseline.py       # classical load flow + solver comparison
```

Each script prints its tables and writes CSVs to `results/tables/`.

**What WP1 should tell you** (confirmed on this machine):

| case | n_bus | NR | GS | FDXB | FDBX |
|---|---|---|---|---|---|
| case9 | 9 | 4 | 74 | 10 | 9 |
| case14 | 14 | 4 | 77 | 10 | 13 |
| case30 | 30 | 4 | 247 | 14 | 9 |
| case39 | 39 | 4 | **diverged** | 12 | 13 |

Newton-Raphson takes 4 iterations regardless of system size — that is quadratic
convergence. Gauss-Seidel scales terribly (74 → 247) and on IEEE 39-bus fails
outright at acceleration 1.6, converging only at 1.4 or below. **That parameter
sensitivity is the point of the comparison**: NR has no such parameter.

Sanity gates that must pass before trusting anything downstream:
- our Ybus vs pandapower's internal Ybus: `< 1e-13`
- our NR/FDLF vs pandapower voltages: `< 1e-8 pu`
- Gauss-Seidel is gated at `1e-6 pu` because it converges on |ΔV| to 1e-8, not on
  power mismatch — holding it to the Newton gate would fail it for being asked a
  different question

---

## Track B — MATLAB / MATPOWER

### B1. Every session

```matlab
cd('<...>\my-portfolio\IEEE\matlab')
setup_paths
```

### B2. WP1 — the MATPOWER twin

```matlab
wp1_baseline_matlab
```

This does three things:
1. Validates `gb_ybus` against MATPOWER's `makeYbus` (must be `< 1e-9`)
2. Runs **both** our hand-written `gb_newton` / `gb_gauss` / `gb_fdlf` **and**
   MATPOWER's own `runpf` with `pf.alg` set to NR / GS / FDXB / FDBX
3. Sweeps the Gauss-Seidel acceleration factor

Writes `results/tables/wp1_*_matlab.csv`.

**Why write our own solvers when `runpf` exists.** Two reasons, and they are the
reason this project is not a `runpf` wrapper. First, `runpf` is a black box: it
will not tell you its Jacobian condition number, its mismatch history, or how
many iterations Gauss-Seidel needed at a given acceleration factor. Second, an
interviewer asking whether you understand load flow is asking about the
derivation, not the API. `gb_newton` builds the polar Jacobian from
`dS/d|V|` and `dS/dθ` explicitly.

### B3. What to check

| Check | Expected |
|---|---|
| `gb_ybus` vs `makeYbus` | `< 1e-9` |
| `gb_newton` vs `runpf` | `< 1e-8 pu` |
| `gb_newton` iterations | 4, on every case |
| `gb_gauss` on case39 @ accel 1.6 | does not converge — this is correct |

---

## Track C — Stateflow (the ride-through logic)

```matlab
setup_paths
build_frt_chart
open_system('gb_frt_chart')
```

Builds the IEEE 2800 fault-ride-through state machine:

```
NORMAL ──[V < 0.88]──> LVRT ──[V > 0.90]──> RECOVERY ──> NORMAL
   │                     │
   └──[V > 1.10]──> HVRT │
                         └──[after(t_env, sec) & V < v_env]──> TRIPPED
```

**Why Stateflow and not a Simulink block diagram.** This behaviour is discrete
and event-driven: mode changes, entry actions, a latching trip, and a timer whose
limit depends on how deep the voltage dip is. Modelling that with switches and
relational operators produces something unreadable and hard to prove correct.
Stateflow is the right tool, and the chart is close to a direct transcription of
the clause in the standard.

The critical behaviour is the **priority flip**. In `NORMAL` the converter runs
active-power priority (you are paid for MWh). On entering `LVRT` it switches to
**reactive priority** and injects `Iq = K1·(1 − V)` — this is what IEEE 2800
requires, and it is the mechanism that keeps negative-sequence current available
for the protection relays downstream.

The ride-through envelope is read from `config/scenarios.yaml` (`ibr.lvrt`), so
the chart and the Python fault solver enforce the same curve.

### Verify it

Feed `Vpu` a step 1.0 → 0.2 → 1.0:
- a **120 ms** dip → `state_id` goes 0 → 1 → 3 → 0, never 4 (inside the 150 ms
  zero-voltage envelope, so it must ride through)
- a **200 ms** dip → `state_id` reaches 4 (TRIPPED)

If a 120 ms dip trips, the envelope is wired wrong.

---

## Track D — Simulink converter models

```matlab
build_converter_control('gfl')     % grid-following
build_converter_control('gfm')     % grid-forming
```

These build the **control** side from plain Simulink blocks. Each subsystem
carries its own specification as an annotation — the maths is in the model.

**Why control and network are built separately.** Control logic scripted this way
is portable and unit-testable on its own. The Simscape electrical network is far
easier and more reliable to wire in the GUI, and scripting it would produce a
brittle script that breaks between releases. So: script the control, draw the
network.

### The GFL / GFM distinction — the thing to actually understand

| | Grid-following (GFL) | Grid-forming (GFM) |
|---|---|---|
| Synchronisation | PLL tracks grid angle | Sets its own angle from power droop |
| Behaves as | Current source | Voltage source behind an impedance |
| Weak grid | Unstable as SCR falls — PLL sees its own current move the terminal voltage | Stable at SCR < 1 |
| Black start | Cannot | Can |
| Inertia | None | Synthetic, if VSM control is enabled |
| Current limiting | Hard clamp | Virtual impedance |

The virtual-impedance detail matters: a GFM converter that simply clamps its
current reverts to grid-following behaviour during the fault and loses the
advantage you built it for.

### Building the Simscape network (GUI, ~30 min for IEEE 9-bus)

1. `powergui` block — **required**, add it first. Set **Discrete, Ts = 50e-6**.
2. Three-Phase Source for the slack, with the Thevenin impedance from the
   MATPOWER Zbus (Python prints it: `ybus.thevenin_impedance`).
3. Three-Phase PI Section Lines from the branch R/X/B.
4. Three-Phase Series RLC Loads, "constant Z" type, from the case bus data.
5. Synchronous Machine pu Standard for machines that stay synchronous;
   Average-Value Inverter + your control model for the ones displaced by IBRs.
6. Three-Phase Fault block at the fault bus, with transition times from
   `config/scenarios.yaml` (`dynamics.fault_apply_s`, `fault_duration_sweep_s`).
7. Three-Phase V-I Measurement blocks at each relay point.

**Validate before trusting any dynamic result:** run the model to steady state
and check the bus voltages against the MATPOWER load flow for the same case. If
the pre-fault operating point does not match to ~1e-3 pu, the network is wired
wrong and every transient result from it is meaningless.

---

## Track E — Simscape network (the EMT model)

```matlab
setup_paths
build_ieee9_simscape            % discrete, 50 us
% ... wire the blocks in the GUI ...
validate_ieee9_steady_state     % GATE -- must pass before anything else
```

`build_ieee9_simscape` parameterises every block from MATPOWER `case9`: line
R/L/C from the branch data, loads from the bus data, source impedances from the
same assumed short-circuit level the phasor fault studies use. **It does not
draw the wires** — connections are far more reliable to make in the GUI than to
script, and a scripted wiring routine breaks between releases.

`validate_ieee9_steady_state` is not optional. A transient simulation answers
"what happens after the disturbance"; if the operating point it starts *from* is
wrong, so is the answer, at any timestep. The gate compares the settled model
against `gb_newton` on the same case and fails with a ranked list of likely
causes.

---

## Track F — Transient stability and CCT

```matlab
setup_paths
wp7_transient_matlab
```

**Read the CCT warning before interpreting anything from this.** Critical
clearing time is *not* comparable across scenarios with different numbers of
dynamic units. The criterion is rotor-angle separation, which only exists for
units that have a rotor angle — so converting machines to grid-following
converters makes CCT look better by deleting the failure mode being measured.
Measured on IEEE 39-bus: 0.121 s at 0% IBR "improves" to 0.977 s at 40.6%
grid-following penetration. That is an artefact.

The controlled comparison is the virtual-inertia sweep (part 3): fixed topology,
fixed dispatch, fixed unit count, only the grid-forming inertia varies. That
gives **4.1× longer CCT** across H = 0.5 → 16 s, and it is where the
grid-forming claim is actually settled.

Baseline check: 0.121 s on the intact IEEE 39-bus against a published ~0.15 s.
The machine data comes from `machine_dynamics` in `scenarios.yaml` (Anderson &
Fouad for case9, Athay et al. for case39). A generic inertia constant instead
gives a CCT roughly 4× too long — this is one place where invented data would
quietly wreck the result.

---

## Track G — the remaining Stateflow charts

```matlab
build_current_limiter_chart     % priority selection under the current limit
build_relay_chart               % pickup / time / trip / reclose / lockout
```

**Current limiter.** When the commanded current exceeds the limit, something has
to give. Active priority keeps P (right in normal operation); reactive priority
keeps Q (required during a fault by IEEE 2800, because it is what stops the
voltage collapsing further). The chart switches mode on the FRT chart's state
output — that is the coupling that makes the pair standards-compliant.

Verify: drive `Id_cmd = 1.0, Iq_cmd = 0.8` (magnitude 1.28 > 1.20).
- ACTIVE_PRIORITY → Id = 1.00, Iq = 0.663
- REACTIVE_PRIORITY → Iq = 0.80, Id = 0.894

**Relay.** The full sequence, including the reclose cycle and lockout count that
a single-shot calculation cannot express. Verify three behaviours:
1. `Imag = 3 pu, Zapp = 0.5` → zone 1, immediate trip, reclose
2. `Imag = 3 pu` sustained → two shots, then LOCKOUT
3. `Imag = 1.1 pu` → **never leaves MONITOR** — this is the WP4 failure mode:
   inverter-limited fault current below pickup means the relay never sees the
   fault at all

---

## Order of work

```
1. Python: export_config.py            (generates the MATLAB bridge too)
2. Python: wp1_baseline.py             (baseline + gates)
3. MATLAB: setup_paths, wp1_baseline_matlab
4. Compare the two WP1 tables          <- if they disagree, stop and fix
5. Python: WP2..WP6 studies
6. MATLAB: WP2/WP3 twins
7. Stateflow: build_frt_chart, verify the 120 ms / 200 ms dip behaviour
8. Simulink: build converter control, wire the Simscape network, validate
   the pre-fault operating point against MATPOWER
9. Dynamic studies: CCT search, GFM vs GFL
```

Never skip step 4. Everything after it assumes the two toolchains agree on the
classical problem.
