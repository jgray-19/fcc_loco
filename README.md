# fcc_loco

Simulated LOCO on FCC-ee (LCC_107, Z), 4D, no radiation, fitted with accelmo's POCO closed-orbit fitter.

- **Truth and targets: xsuite.** Random relative k1L errors (0.5e-2 rms) on every arc quadrupole, tunes re-matched. Targets are the BPM orbit changes for +-0.5 urad trims of N H and N V correctors spread around the ring, plus 0.1 um BPM noise.
- **Fit: MAD-NG.** `accelmo.poco.ClosedOrbitFitter` with one `dk1l` knob per arc quadrupole (`accelmo.machine.accelerators.fcc.FCC`).
- **Output:** `beta_beat.png` (before and after, against the nominal xsuite line), `dk1l_correlation.png` (fitted against actual), `metrics.json`.

Layout: `fcc_loco/lattice.py` (xsuite line, errors, corrector choice), `scan.py` (targets), `fit.py` (the accelmo call, about 30 lines), `analysis.py`, `plots.py`; `scripts/run_loco.py` runs it all.

## Run

Heavy (one MAD-NG process per trim, 2639 BPMs): run on a compute host, not a login node.

```
scripts/00_fetch_lattice.sh
uv run scripts/00_make_model.py
uv run scripts/run_loco.py --n-correctors 10 --seed 1
```

Results go to `results/loco_n<N>_seed<S>/`. `uv run pytest` runs the xsuite-only checks.
