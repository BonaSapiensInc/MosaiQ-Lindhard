#!/usr/bin/env python3
"""Render Static Structure Factors S(q) for various Gamma values (PRE style)."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILE = PROJECT_ROOT / "output" / "output_Sq_gamma.dat"
OUTPUT_PDF = PROJECT_ROOT / "output" / "Sq_gamma_sweep.pdf"
OUTPUT_PNG = PROJECT_ROOT / "output" / "Sq_gamma_sweep.png"

COLOR_EE = "#1f3b73"
COLOR_II = "#8b1e1e"


def configure_seaborn() -> None:
    sns.set_theme(
        style="whitegrid",
        context="paper",
        font_scale=1.1,
        rc={
            "font.family": "serif",
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "axes.linewidth": 0.8,
            "grid.alpha": 0.3,
            "savefig.dpi": 300,
        },
    )


def parse_gamma_data(path: Path) -> dict[float, dict[str, np.ndarray]]:
    """Parse the output_Sq_gamma.dat file containing blocks of Gamma sweeps."""
    if not path.is_file():
        raise SystemExit(f"Error: data file not found: {path}")

    data_by_gamma: dict[float, dict[str, list[float]]] = defaultdict(
        lambda: {"q": [], "S_ee": [], "S_ii": []}
    )
    current_gamma: float | None = None

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            gamma_match = re.search(r"Gamma\s*=\s*([0-9.]+)", line)
            if gamma_match:
                current_gamma = float(gamma_match.group(1))
                continue

            if line.startswith("#"):
                continue

            if current_gamma is not None:
                parts = line.split()
                if len(parts) >= 3:
                    data_by_gamma[current_gamma]["q"].append(float(parts[0]))
                    data_by_gamma[current_gamma]["S_ee"].append(float(parts[1]))
                    data_by_gamma[current_gamma]["S_ii"].append(float(parts[2]))

    return {
        gamma: {
            "q": np.array(block["q"]),
            "S_ee": np.array(block["S_ee"]),
            "S_ii": np.array(block["S_ii"]),
        }
        for gamma, block in data_by_gamma.items()
    }


def render_plot(data: dict[float, dict[str, np.ndarray]]) -> None:
    configure_seaborn()

    gammas = sorted(data.keys())
    if len(gammas) != 4:
        print(f"Warning: Expected 4 Gamma values for a 2x2 grid, found {len(gammas)}")

    fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for idx, gamma in enumerate(gammas):
        if idx >= 4:
            break
        ax = axes_flat[idx]
        q = data[gamma]["q"]
        s_ee = data[gamma]["S_ee"]
        s_ii = data[gamma]["S_ii"]

        ax.plot(q, s_ee, color=COLOR_EE, linewidth=2.0, label=r"$S_{ee}(q)$")
        ax.plot(q, s_ii, color=COLOR_II, linewidth=2.0, label=r"$S_{ii}(q)$")

        ax.set_title(rf"$\Gamma = {gamma:g}$", fontweight="bold")
        ymax = max(float(np.max(s_ee)), float(np.max(s_ii)), 0.1)
        ax.set_ylim(0.0, ymax * 1.1)
        ax.set_xlim(0, float(np.max(q)))

        if idx in (2, 3):
            ax.set_xlabel(r"Reduced wavevector $\bar{q} = q/k_F$")
        if idx in (0, 2):
            ax.set_ylabel(r"Static Structure Factor $S(\bar{q})$")
        if idx == 0:
            ax.legend(loc="upper left", frameon=True, edgecolor="0.7")

    plt.tight_layout()
    OUTPUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PDF)
    fig.savefig(OUTPUT_PNG)
    plt.close(fig)
    print(f"Saved {OUTPUT_PDF} and {OUTPUT_PNG}")


def main() -> None:
    data = parse_gamma_data(DATA_FILE)
    render_plot(data)


if __name__ == "__main__":
    main()
