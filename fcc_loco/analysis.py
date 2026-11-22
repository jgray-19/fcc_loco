"""Score a fit against the truth: knob errors and beta beating at the BPMs (xsuite)."""

import numpy as np
import pandas as pd

from fcc_loco.lattice import apply_dk1l, load_line, twiss_at_bpms


def knob_errors(fit: dict[str, float], truth: dict[str, float]) -> dict[str, float]:
    names = sorted(truth)
    t = np.array([truth[n] for n in names])
    f = np.array([fit.get(n, 0.0) for n in names])
    return {
        "dk1l_truth_rms": float(np.sqrt(np.mean(t**2))),
        "dk1l_error_rms": float(np.sqrt(np.mean((f - t) ** 2))),
        "dk1l_correlation": float(np.corrcoef(f, t)[0, 1]),
        "dk1l_slope": float(np.polyfit(t, f, 1)[0]),
    }


def _beat(twiss: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "s": reference["s"],
            "bbx": twiss["betx"] / reference["betx"] - 1,
            "bby": twiss["bety"] / reference["bety"] - 1,
        }
    )


def beta_beats(truth: dict[str, float], fitted: dict[str, float]):
    """``(before, after)``: beta beating against the nominal line of the perturbed line and of the
    perturbed line with the fitted errors removed."""
    nominal = twiss_at_bpms(load_line())
    before = twiss_at_bpms(apply_dk1l(load_line(), truth))
    residual = {name: error - fitted.get(name, 0.0) for name, error in truth.items()}
    after = twiss_at_bpms(apply_dk1l(load_line(), residual))
    return _beat(before, nominal), _beat(after, nominal)


def beat_rms(table: pd.DataFrame) -> dict[str, float]:
    return {c: float(np.sqrt(np.mean(table[c].to_numpy() ** 2))) for c in ("bbx", "bby")}
