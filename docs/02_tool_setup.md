# Tool Setup

Everything you need to install, once, before running anything. Ordered by how
much you need it: the Python stack runs every study, MATLAB adds the
independent cross-check and the dynamic work.

Verified working on this machine: **Python 3.13.7**, **MATLAB R2025b**, Windows 11.

---

## 1. Python + pandapower (required)

Runs every phasor-domain study in this project. No licence, no cost.

```bash
cd IEEE/python
python -m pip install -r ../requirements.txt
```

Verify:

```bash
python -c "import pandapower as pp; print(pp.__version__)"     # expect >= 3.5
python studies/wp1_baseline.py                                  # should print the WP1 tables
```

**Installed and confirmed working here:** pandapower 3.5.4, numpy 2.4.2,
numba 0.66.0 (numba is optional but pandapower is noticeably slower without it).

> The IEEE test cases ship inside pandapower (`pandapower.networks.case9/14/30/39`),
> so there is nothing to download. `gridbench.ppc.load_ppc()` pulls them out in
> MATPOWER format, which is what makes the Python and MATLAB sides comparable.

---

## 2. MATLAB + MATPOWER (required for the cross-check)

MATPOWER is a free add-on; it is not bundled with MATLAB.

1. Download MATPOWER 8.x from <https://matpower.org>
2. Unzip somewhere permanent, e.g. `C:\matpower8`
3. In MATLAB:
   ```matlab
   addpath(genpath('C:\matpower8'));
   savepath
   ```
4. Verify: `test_matpower` — expect all tests to pass (takes a few minutes)

Then, in every MATLAB session for this project:

```matlab
cd('<...>\my-portfolio\IEEE\matlab')
setup_paths          % adds subfolders, checks MATPOWER, checks the config bridge
```

`setup_paths` tells you exactly what is missing if anything is.

### The config bridge (do this before your first MATLAB run)

MATLAB and Python read the **same** `config/scenarios.yaml`. MATLAB's YAML
support varies by release, so the YAML is mechanically translated to JSON:

```bash
cd IEEE/python
python studies/export_config.py
```

Re-run it any time you edit `config/*.yaml`. If you forget, `gb_config` raises an
error telling you to — it will not silently use stale values.

---

## 3. Simulink + Simscape Electrical + Stateflow (for the dynamic work)

Needed only for WP9–WP12 (converter models, ride-through logic, CCT search).
Check what you actually have:

```matlab
ver          % look for Simulink, Simscape Electrical, Stateflow
```

Or test each specifically:

```matlab
license('test','Simulink')
license('test','Simscape_Electrical')
license('test','Stateflow')
```

If Simscape Electrical is missing, everything in WP1–WP8 still runs — the
phasor-domain work does not depend on it.

---

## 4. What is NOT set up yet (deliberately deferred)

| Tool | Status | Constraint to plan around |
|---|---|---|
| PSCAD | deferred | Free Edition capped at **15 nodes** → needs a Thevenin-reduced POI equivalent |
| ETAP | deferred | Student Edition capped at **15 buses** → IEEE 9-bus only |
| DIgSILENT PowerFactory | deferred | PF4S student licence: **50 nodes**, ~€150/year |
| OpenDSS | deferred | Free (`pip install py-dss-interface`); for the IEEE 13/123-node feeders |

The `config/scenarios.yaml` entries for these already exist so the work can be
picked up without restructuring anything.

---

## Version pins

Recorded so results stay reproducible. If a number in the report cannot be
regenerated, the first thing to check is whether these moved.

| Component | Version used |
|---|---|
| Python | 3.13.7 |
| pandapower | 3.5.4 |
| numpy | 2.4.2 |
| numba | 0.66.0 |
| MATLAB | R2025b |
| MATPOWER | 8.x |
