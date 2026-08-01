# Suggested portfolio entry

Your `src/data/content.js` currently lists this project as **planned**, with the
textbook description. Below is a replacement entry built from the actual results.

**Status: `content.js` HAS now been updated** (project entry + skills), with your
go-ahead. The variants below are kept for reference if you want to reword.

---

## Replacement for the entry at `src/data/content.js:210-220`

```js
{
  title: "IEEE Bus System — IBR Penetration Benchmark (Load Flow, Faults, Stability & Security)",
  status: "published",
  note: "Cross-tool power systems benchmark",
  description:
    "Quantifies where classical power system analysis stops being valid as inverters " +
    "displace synchronous machines, on IEEE 9/14/30/39-bus systems. Two independent " +
    "toolchains (Python+pandapower, MATLAB+MATPOWER) solve identical problems and are " +
    "diffed against each other, with Simulink/Simscape/Stateflow for converter control. " +
    "Key results: IEEE 39-bus breaches the SCR≥3 interconnection screen at 28.9% IBR " +
    "penetration; fault current falls 62% and relay misoperation reaches 91% at 84% " +
    "penetration; RoCoF grows from 1.0 to 6.3 Hz/s. Demonstrates two mitigations — " +
    "25% grid-forming control restores a solvable fault problem at 100% inverter " +
    "penetration and 4.1× critical clearing time, and synchronous condensers lift " +
    "short-circuit strength 162%. Includes an IEEE 2800-2022 ride-through Stateflow " +
    "implementation and a demonstration that exact AC false-data injection is " +
    "structurally undetectable by chi-squared residual testing. 56 regression gates, " +
    "22 result tables, 30 figures, fully reproducible from one command.",
  tags: [
    "MATPOWER", "pandapower", "Simulink", "Simscape", "Stateflow",
    "IEEE 2800", "Grid-Forming", "Short-Circuit Ratio", "Transient Stability",
    "Protection", "State Estimation", "Python", "MATLAB",
  ],
  githubUrl: "https://github.com/MukeshSingh123-tech/ieee-ibr-benchmark",
  liveUrl: "",
},
```

## A shorter variant, if the card is length-constrained

```js
{
  title: "IEEE Bus System — IBR Penetration Benchmark",
  status: "published",
  note: "Cross-tool power systems benchmark",
  description:
    "Measures where classical load flow, symmetrical-component fault analysis and " +
    "conventional protection stop being trustworthy as inverters displace synchronous " +
    "machines on IEEE 9/14/30/39-bus systems. Python+pandapower and MATLAB+MATPOWER " +
    "cross-validated to 1e-8 pu, plus Simulink/Stateflow converter control. " +
    "IEEE 39-bus fails the SCR≥3 screen at 28.9% penetration; relay misoperation hits " +
    "91% at 84%. Shows grid-forming control and synchronous condensers both fix it, " +
    "with a 4.1× gain in critical clearing time. 56 regression gates, reproducible " +
    "end-to-end.",
  tags: ["MATPOWER", "pandapower", "Simulink", "Stateflow", "IEEE 2800",
         "Grid-Forming", "Protection", "Python", "MATLAB"],
  githubUrl: "https://github.com/MukeshSingh123-tech/ieee-ibr-benchmark",
  liveUrl: "",
},
```

---

## Before you publish

1. **Create the GitHub repo** — `content.js` now points at
   `github.com/MukeshSingh123-tech/ieee-ibr-benchmark`, which does NOT exist yet
   and will 404 until you create and push it. The local repo is committed and
   the remote is already configured.
2. **Consider what to commit.** `results/tables/` and `results/figures/` are
   small and are the evidence for every claim, so they are worth committing even
   though they are generated.
3. **The README is the landing page.** It leads with results and states
   limitations explicitly — that combination is what makes it read as research
   rather than a coursework write-up.

## Novelty check — corrected

An earlier draft of this file claimed two findings were novel. **A literature
check does not support that.** Recorded here so the claim is not repeated:

| Finding | Status | Evidence |
|---|---|---|
| Rotor-angle CCT is invalid for IBR-penetration comparisons | **already published** | IET, *Application of an Advanced Short Circuit Strength Metric to Evaluate Ireland's High Renewable Penetration Scenarios* — states that metrics including CCT "may not be appropriate as metrics of system strength in systems with high levels of IBRs because of the absence of rotating inertia" |
| Classical fixed-Qmax error flips sign with operating point | **mechanism already known** | WECC *IBR Power Plant Modeling and Validation Guideline* — actual Qmax/Qmin "is limited by the power factor at the active power output"; the P-Q capability being a curve, not a constant, is established practice |
| Exact AC FDIA defeats chi-squared detection | **reproduces Liu, Ning & Reiter (2009)** | the foundational FDIA paper |

**So do not claim novelty.** What this project actually offers, which is still
worth presenting:

- **Quantification on standard benchmarks** — most of the above is stated
  qualitatively in the literature. Here it is measured: SCR breach at 28.9%,
  fault current −62%, misoperation 91%, RoCoF 1.0 → 6.3 Hz/s, CCT gain 4.1×.
- **Reproducibility** — one command regenerates every table and figure, 56
  regression gates, hypotheses recorded before running with outcomes written
  back after.
- **A controlled experiment design** — the inverter is sized so its capability
  exactly matches the machine it replaces (verified to 1e-14), so the result is
  attributable to the physics and not to a sizing assumption. That discipline is
  not universal in the literature.
- **Two independent toolchains cross-validated to 1e-8 pu** on the same problems.
- **The AC-exact vs DC-linearised FDIA contrast** — a cleaner demonstration than
  the original, and it connects directly to the ChainPMU work.

Frame it as *rigorous engineering and a reproducible benchmark*, not as new
science. That is a claim you can defend in front of a power systems engineer.
