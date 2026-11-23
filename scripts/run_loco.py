# /// script
# requires-python = ">=3.11"
# dependencies = ["fcc-loco"]
# [tool.uv.sources]
# fcc-loco = { path = "..", editable = true }
# ///
"""FCC-ee LOCO: xsuite perturbed line gives the targets, MAD-NG fits one k1 knob per arc quadrupole.

Writes beta_beat.png, dk1l_correlation.png and metrics.json to results/loco_<tag>/.
"""

import argparse
import json
import os

import numpy as np

from fcc_loco.analysis import beat_rms, beta_beats, knob_errors
from fcc_loco.config import BPM_NOISE, RESULTS_DIR
from fcc_loco.fit import build_series, run_fit
from fcc_loco.lattice import draw_truth, load_line, select_correctors
from fcc_loco.plots import plot_beta_beat, plot_knob_correlation
from fcc_loco.scan import measured, truth_delta_orbits

parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=int, default=1)
parser.add_argument(
    "--n-correctors", type=int, default=10, help="per plane, spread around the ring"
)
parser.add_argument("--noise", type=float, default=BPM_NOISE, help="rms BPM noise [m]")
parser.add_argument(
    "--jobs", type=int, default=min(8, os.cpu_count()), help="xsuite scan processes"
)


def main() -> None:
    args = parser.parse_args()

    out = RESULTS_DIR / f"loco_n{args.n_correctors}_seed{args.seed}"
    out.mkdir(parents=True, exist_ok=True)

    correctors = select_correctors(load_line(), args.n_correctors)
    truth = draw_truth(args.seed)
    delta = truth_delta_orbits(truth, correctors, args.jobs)
    targets = measured(delta, np.random.default_rng(3000 + args.seed), args.noise)

    fitted, losses = run_fit(
        build_series(correctors, targets),
        correctors,
        out,
    )

    before, after = beta_beats(truth, fitted)
    plot_beta_beat(before, after, out / "beta_beat.png")
    plot_knob_correlation(truth, fitted, out / "dk1l_correlation.png")
    metrics = {
        "iterations": len(losses),
        "beta_beat_before": beat_rms(before),
        "beta_beat_after": beat_rms(after),
        "knobs": knob_errors(fitted, truth),
    }
    (out / "metrics.json").write_text(json.dumps(metrics, indent=1))
    print(json.dumps(metrics, indent=1))


if __name__ == "__main__":
    main()
