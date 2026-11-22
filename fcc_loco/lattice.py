"""Xsuite FCC-ee line, random quadrupole errors and the corrector selection."""

import json
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
import xtrack as xt

from fcc_loco.config import ARC_QUAD, QUAD_REL_DK1, XSUITE_JSON, XSUITE_LINE

CORRECTOR = r"^[hv]cor_"


@dataclass(frozen=True)
class Corrector:
    name: str
    plane: str  # "x" or "y"
    s: float  # m

    @property
    def variable(self) -> str:
        return f"k_{self.name}"


def load_line() -> xt.Line:
    """The line with cavities off and full bend edges (as MAD-NG's dipole fringe)."""
    dct = json.loads(XSUITE_JSON.read_text())
    valid = set(xt.Quadrupole.get_available_integrators())
    for element in dct["elements"].values():
        if "integrator" in element and element["integrator"] not in valid:
            element["integrator"] = "yoshida4"
    line = xt.Environment.from_dict(dct).lines[XSUITE_LINE]
    for element in line.element_dict.values():
        if isinstance(element, xt.Cavity):
            element.voltage = 0.0
        elif isinstance(element, xt.RBend | xt.Bend):
            element.edge_entry_model = "full"
            element.edge_exit_model = "full"
    return line


def arc_quadrupoles(line: xt.Line) -> dict[str, tuple[float, float]]:
    """``{name: (k1, length)}`` of the arc quadrupoles."""
    return {
        name: (float(element.k1), float(element.length))
        for name in line.element_names
        if re.match(ARC_QUAD, name)
        and type(element := line.element_dict[name]).__name__ == "Quadrupole"
        and element.k1 != 0.0
    }


def apply_dk1l(line: xt.Line, dk1l: dict[str, float]) -> xt.Line:
    """Add ``dk1l`` (keys ``<quad>.dk1l``) to the quadrupoles' integrated strengths, in place."""
    for key, value in dk1l.items():
        element = line.element_dict[key.rsplit(".", 1)[0]]
        element.k1 = float(element.k1) + value / float(element.length)
    return line


def _match_tunes(line: xt.Line, qx: float, qy: float, quads: list[str]):
    """Restore the tunes of the perturbed ``line`` with one QF and one QD trim knob; returns the
    trims as ``dk1l``."""
    base = {
        name: (float(line.element_dict[name].k1), float(line.element_dict[name].length))
        for name in quads
    }
    for knob in ("dkqf_trim", "dkqd_trim"):
        line.vars[knob] = 0.0
    for name, (k1, _) in base.items():
        knob = "dkqf_trim" if name.startswith("qf") else "dkqd_trim"
        line.element_refs[name].k1 = k1 * (1.0 + line.vars[knob])
    line.match(
        method="4d",
        vary=xt.VaryList(["dkqf_trim", "dkqd_trim"], step=1e-6),
        targets=xt.TargetSet(qx=qx, qy=qy, tol=1e-7),
    )
    trims = {"f": float(line.vars["dkqf_trim"]._value), "d": float(line.vars["dkqd_trim"]._value)}
    return {
        f"{name}.dk1l": k1 * trims["f" if name.startswith("qf") else "d"] * length
        for name, (k1, length) in base.items()
    }


def draw_truth(seed: int, rel_error: float = QUAD_REL_DK1) -> dict[str, float]:
    """``{<quad>.dk1l: error}``, sigma ``rel_error * |k1L|``, plus the tune-restoring trims."""
    line = load_line()
    quads = arc_quadrupoles(line)
    rng = np.random.default_rng(seed)
    errors = {
        f"{name}.dk1l": float(rng.normal(0, rel_error * abs(k1 * length)))
        for name, (k1, length) in quads.items()
    }
    nominal = line.twiss(method="4d")
    trims = _match_tunes(
        apply_dk1l(load_line(), errors), float(nominal.qx), float(nominal.qy), list(quads)
    )
    return {key: value + trims[key] for key, value in errors.items()}


def select_correctors(line: xt.Line, n_per_plane: int) -> list[Corrector]:
    """``n_per_plane`` H and V correctors at evenly spaced s, next to arc quadrupoles."""
    candidates = [
        name
        for name in line.element_names
        if re.match(CORRECTOR, name) and re.match(ARC_QUAD, name.split("_", 1)[1])
    ]
    table = line.get_table()
    circumference = float(line.get_length())
    chosen = []
    for prefix, plane in (("h", "x"), ("v", "y")):
        pool = [Corrector(n, plane, float(table["s", n])) for n in candidates if n[0] == prefix]
        positions = np.array([c.s for c in pool])
        for target in (np.arange(n_per_plane) + 0.5) * circumference / n_per_plane:
            chosen.append(pool[int(np.argmin(np.abs(positions - target)))])
    return chosen


def set_kick(line: xt.Line, corrector: Corrector, kick: float) -> None:
    """MAD-X kick convention: ``hkick`` is ``knl[0] = -kick``, ``vkick`` is ``ksl[0] = kick``."""
    element = line.element_dict[corrector.name]
    if corrector.plane == "x":
        element.knl[0] = -kick
    else:
        element.ksl[0] = kick


def twiss_at_bpms(line: xt.Line) -> pd.DataFrame:
    """x, y, betx, bety, s at every BPM (4D, on momentum), index = BPM name."""
    names = [name for name in line.element_names if name.startswith("bpm_")]
    rows = line.twiss(method="4d", delta0=0.0, chrom=False).rows[names]
    columns = ("x", "y", "betx", "bety", "s")
    return pd.DataFrame({c: np.asarray(rows[c]) for c in columns}, index=list(rows.name))
