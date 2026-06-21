#!/usr/bin/env python3
"""Publication-quality plasmon dispersion and Landau damping figure (PRE style)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILE = PROJECT_ROOT / "output" / "output_plasmon_dispersion.dat"
OUTPUT_PDF = PROJECT_ROOT / "output" / "plasmon_dispersion.pdf"
OUTPUT_PNG = PROJECT_ROOT / "output" / "plasmon_dispersion.png"

COLOR_OMEGA_P = "#1f3b73"  # indigo
COLOR_BOHM_GROSS = "#4d4d4d"
COLOR_LANDAU = "#8b1e1e"  # crimson


def load_dispersion(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if not path.is_file():
        raise SystemExit(f"Error: data file not found: {path}")

    table = np.loadtxt(path, comments="#")
    if table.ndim != 2 or table.shape[1] < 4:
        raise SystemExit(f"Error: expected 4 columns in {path}, got shape {table.shape}")

    q = table[:, 0]
    omega_p = table[:, 1]
    landau = table[:, 2]
    bohm_gross = table[:, 3]
    return q, omega_p, landau, bohm_gross


def finite_segments(q: np.ndarray, y: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split (q, y) into contiguous segments where y is finite."""
    mask = np.isfinite(y)
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    start: int | None = None

    for index, is_finite in enumerate(mask):
        if is_finite and start is None:
            start = index
        elif not is_finite and start is not None:
            segments.append((q[start:index], y[start:index]))
            start = None

    if start is not None:
        segments.append((q[start:], y[start:]))

    return segments


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 0.8,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.5,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_dispersion(
    q: np.ndarray,
    omega_p: np.ndarray,
    landau: np.ndarray,
    bohm_gross: np.ndarray,
) -> plt.Figure:
    configure_matplotlib()

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(6.5, 5.5),
        gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.08},
    )

    for q_seg, y_seg in finite_segments(q, omega_p):
        ax_top.plot(
            q_seg,
            y_seg,
            color=COLOR_OMEGA_P,
            linewidth=1.8,
            solid_capstyle="round",
        )

    ax_top.plot(
        q,
        bohm_gross,
        color=COLOR_BOHM_GROSS,
        linewidth=1.4,
        linestyle="--",
        dashes=(5, 3),
        label="Bohm--Gross",
    )

    ax_top.set_ylabel(r"$\bar{\omega}_p$ [$E_\mathrm{F}/\hbar$]")
    ax_top.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.35)
    ax_top.legend(
        [
            plt.Line2D([0], [0], color=COLOR_OMEGA_P, linewidth=1.8),
            plt.Line2D([0], [0], color=COLOR_BOHM_GROSS, linewidth=1.4, linestyle="--"),
        ],
        [r"Extracted $\bar{\omega}_p$", "Bohm--Gross"],
        loc="upper left",
        frameon=True,
        framealpha=0.92,
        edgecolor="0.75",
    )

    for q_seg, y_seg in finite_segments(q, landau):
        positive = y_seg > 0.0
        if np.count_nonzero(positive) >= 2:
            ax_bottom.plot(
                q_seg[positive],
                y_seg[positive],
                color=COLOR_LANDAU,
                linewidth=1.8,
                solid_capstyle="round",
            )

    ax_bottom.set_yscale("log")
    ax_bottom.set_ylabel(r"$\Im[\varepsilon(\bar{q}, \bar{\omega}_p)]$")
    ax_bottom.set_xlabel(r"Reduced wavevector $\bar{q} = q/k_\mathrm{F}$")
    ax_bottom.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.35)

    fig.align_ylabels([ax_top, ax_bottom])
    fig.subplots_adjust(hspace=0.08)
    return fig


def main() -> None:
    q, omega_p, landau, bohm_gross = load_dispersion(DATA_FILE)
    fig = plot_dispersion(q, omega_p, landau, bohm_gross)

    OUTPUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PDF)
    fig.savefig(OUTPUT_PNG)
    plt.close(fig)

    n_roots = int(np.sum(np.isfinite(omega_p)))
    print(f"Loaded {q.size} q points ({n_roots} finite plasmon roots).")
    print(f"Saved {OUTPUT_PDF}")
    print(f"Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
