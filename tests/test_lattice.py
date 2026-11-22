"""Physical checks of the xsuite side of the example (no MAD-NG)."""

import re

import numpy as np
import pytest

from fcc_loco.config import ARC_QUAD
from fcc_loco.lattice import (
    apply_dk1l,
    arc_quadrupoles,
    draw_truth,
    load_line,
    select_correctors,
    twiss_at_bpms,
)

QUAD_REL_RMS = 0.5e-2
TUNE_TOL = 1e-6


@pytest.fixture(scope="module")
def line():
    return load_line()


def test_truth_restores_nominal_tunes(line):
    nominal = line.twiss(method="4d")
    perturbed = apply_dk1l(load_line(), draw_truth(1))
    tw = perturbed.twiss(method="4d")
    assert abs(tw.qx - nominal.qx) < TUNE_TOL
    assert abs(tw.qy - nominal.qy) < TUNE_TOL


def test_truth_is_a_beta_beat(line):
    nominal = twiss_at_bpms(line)
    truth = draw_truth(1)
    perturbed = twiss_at_bpms(apply_dk1l(load_line(), truth))
    beat = perturbed["betx"] / nominal["betx"] - 1
    assert 1e-3 < np.sqrt(np.mean(beat**2)) < 0.5
    quads = arc_quadrupoles(line)
    relative = [truth[f"{n}.dk1l"] / (k1 * length) for n, (k1, length) in quads.items()]
    assert np.std(relative) == pytest.approx(QUAD_REL_RMS, rel=0.05)


def test_correctors_are_spread_and_next_to_arc_quads(line):
    chosen = select_correctors(line, 5)
    assert len(chosen) == 10
    assert len({c.name for c in chosen}) == 10
    assert all(re.match(ARC_QUAD, c.name.split("_", 1)[1]) for c in chosen)
    circumference = float(line.get_length())
    for plane in "xy":
        s = np.sort([c.s for c in chosen if c.plane == plane])
        assert np.all(np.diff(s) < 0.4 * circumference)
