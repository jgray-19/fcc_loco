"""Fit the arc-quadrupole k1 knobs to the scan with accelmo's closed-orbit fitter."""

from pathlib import Path

import pandas as pd
from accelmo.fitting.config import OutputConfig, SequenceConfig
from accelmo.machine.accelerators.fcc import FCC
from accelmo.optimisers.levenberg_marquardt import LevenbergMarquardtConfig
from accelmo.poco import ClosedOrbitFitter, ClosedOrbitMeasurement, ClosedOrbitSeries

from fcc_loco.config import SEQUENCE_FILE
from fcc_loco.lattice import Corrector

MAX_ITERATIONS = 100
INITIAL_LAMBDA = 1e-2  # 0 lets the first step overshoot and MAD-NG loses the closed orbit
LAMBDA_DECREASE = 0.05  # per accepted step


def build_series(
    correctors: list[Corrector], targets: dict[tuple[str, float], pd.DataFrame]
) -> list[ClosedOrbitSeries]:
    """One series per trimmed corrector. The orbit of the untrimmed model is zero for k1 errors, so
    the targets (orbit differences) are fitted as absolute orbits."""
    by_name = {c.name: c for c in correctors}
    return [
        ClosedOrbitSeries(
            measurements=(ClosedOrbitMeasurement(orbit=orbit),),
            machine_state={by_name[name].variable: step},
            absolute_planes=("x", "y"),
            label=f"{name}{step:+.1e}",
        )
        for (name, step), orbit in sorted(targets.items())
    ]


def run_fit(
    series: list[ClosedOrbitSeries],
    correctors: list[Corrector],
    output: Path,
    loss_tolerance: float,
) -> tuple[dict[str, float], list[float]]:
    """Returns ``{<quad>.dk1l: fitted}`` (lower case) and the loss of every accepted iteration."""
    accelerator = FCC(sequence_file=str(SEQUENCE_FILE), mode="z", errors={"quad": {"k1"}})
    fitter = ClosedOrbitFitter(
        accelerator=accelerator,
        sequence_config=SequenceConfig(magnet_range="$start/$end"),
        series=series,
        # 4D like the truth (cavities off); every corrector at zero kick
        machine_state={**{c.variable: 0.0 for c in correctors}, "rf_on": 0.0},
        lm_config=LevenbergMarquardtConfig(
            max_iterations=MAX_ITERATIONS,
            initial_lambda=INITIAL_LAMBDA,
            lambda_decrease=LAMBDA_DECREASE,
            loss_relative_tolerance=loss_tolerance,
            gradient_converged_value=0.0,
        ),
        output_config=OutputConfig(run_log_root=output / "runlog", include_uncertainty=False),
    )
    try:
        knobs = fitter.run().knobs
    finally:
        fitter.close()
    return {k.lower(): v for k, v in knobs.items()}, [loss for _, loss in fitter.history]
