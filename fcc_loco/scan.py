"""Synthetic LOCO scan: xsuite orbit changes of the perturbed line for +-TRIM_STEP trims."""

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from fcc_loco.config import TRIM_STEP
from fcc_loco.lattice import Corrector, apply_dk1l, load_line, set_kick, twiss_at_bpms


def _solve_chunk(args) -> dict[tuple[str, float], pd.DataFrame]:
    dk1l, correctors = args
    line = apply_dk1l(load_line(), dk1l)
    reference = twiss_at_bpms(line)[["x", "y"]]
    out = {}
    for corrector in correctors:
        for step in (TRIM_STEP, -TRIM_STEP):
            set_kick(line, corrector, step)
            out[(corrector.name, step)] = twiss_at_bpms(line)[["x", "y"]] - reference
        set_kick(line, corrector, 0.0)
    return out


def truth_delta_orbits(
    dk1l: dict[str, float], correctors: list[Corrector], jobs: int
) -> dict[tuple[str, float], pd.DataFrame]:
    """Noise-free ``orbit(trim) - orbit(untrimmed)`` at the BPMs for every (corrector, step)."""
    jobs = max(1, min(jobs, len(correctors)))
    chunks = [correctors[i::jobs] for i in range(jobs)]
    with ProcessPoolExecutor(jobs, mp_context=mp.get_context("spawn")) as pool:
        parts = pool.map(_solve_chunk, [(dk1l, chunk) for chunk in chunks])
    return {key: frame for part in parts for key, frame in part.items()}


def measured(
    delta_orbits: dict[tuple[str, float], pd.DataFrame], rng: np.random.Generator, noise: float
) -> dict[tuple[str, float], pd.DataFrame]:
    """Targets in the fitter's ``X/ERRX/Y/ERRY`` frame; a difference carries both orbits' noise."""
    sigma = np.sqrt(2.0) * noise
    out = {}
    for key, delta in delta_orbits.items():
        n = len(delta)
        out[key] = pd.DataFrame(
            {
                "X": delta["x"].to_numpy() + rng.normal(0, sigma, n),
                "ERRX": np.full(n, sigma),
                "Y": delta["y"].to_numpy() + rng.normal(0, sigma, n),
                "ERRY": np.full(n, sigma),
            },
            index=pd.Index(delta.index.astype(str)).str.upper(),  # MAD-NG's BPM names
        )
    return out
