"""Beta beating before and after the fit, and fitted against actual quadrupole errors."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BEFORE, AFTER = "#2a78d6", "#eb6834"


def _rms(a) -> float:
    return float(np.sqrt(np.mean(np.square(a))))


def plot_beta_beat(before: pd.DataFrame, after: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True, constrained_layout=True)
    for ax, plane, label in ((axes[0], "bbx", "Horizontal"), (axes[1], "bby", "Vertical")):
        for table, colour, name in ((before, BEFORE, "before"), (after, AFTER, "after")):
            ax.plot(
                table["s"] / 1e3,
                100 * table[plane],
                color=colour,
                lw=0.8,
                label=f"{name}, rms {100 * _rms(table[plane]):.2f}%",
            )
        ax.set_ylabel(f"{label} Δβ/β [%]")
        ax.legend(frameon=False, loc="upper right", ncols=2)
        ax.grid(True, lw=0.5)
    axes[1].set_xlabel("s [km]")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_knob_correlation(truth: dict[str, float], fit: dict[str, float], path: Path) -> None:
    names = sorted(truth)
    actual = np.array([truth[n] for n in names])
    estimated = np.array([fit.get(n, 0.0) for n in names])
    fig, ax = plt.subplots(figsize=(5.5, 5.5), constrained_layout=True)
    ax.scatter(actual, estimated, s=6, color=BEFORE, alpha=0.6, linewidths=0)
    lim = 1.05 * max(np.abs(actual).max(), np.abs(estimated).max())
    ax.plot([-lim, lim], [-lim, lim], color="grey", lw=1, ls="--")
    ax.set_aspect("equal")
    ax.set_xlabel("actual dk1l [1/m]")
    ax.set_ylabel("estimated dk1l [1/m]")
    ax.set_title(
        f"r = {np.corrcoef(estimated, actual)[0, 1]:.4f}, "
        f"residual rms {_rms(estimated - actual):.2e} (true {_rms(actual):.2e})",
        fontsize=10,
    )
    fig.savefig(path, dpi=150)
    plt.close(fig)
