#!/usr/bin/env python3
"""Render Static Structure Factors S(q) per channel with Gamma sweeps (PRE style)."""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from plot_common import OUTPUT_DIR, PDF_RCPARAMS, save_figure

DATA_FILE = OUTPUT_DIR / "output_Sq_gamma.dat"


def configure_seaborn() -> None:
    sns.set_theme(
        style="whitegrid",
        context="paper",
        font_scale=1.1,
        rc={
            "font.family": "serif",
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "axes.linewidth": 0.8,
            "grid.alpha": 0.3,
            **PDF_RCPARAMS,
        },
    )


def parse_gamma_data(path: Path) -> dict[float, dict[str, np.ndarray]]:
    if not path.is_file():
        raise SystemExit(f"Error: data file not found: {path}")

    data_by_gamma: dict[float, dict[str, list[float]]] = defaultdict(
        lambda: {"q": [], "S_ee": [], "S_ii": [], "S_ei": []}
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
                if len(parts) >= 4:
                    data_by_gamma[current_gamma]["q"].append(float(parts[0]))
                    data_by_gamma[current_gamma]["S_ee"].append(float(parts[1]))
                    data_by_gamma[current_gamma]["S_ii"].append(float(parts[2]))
                    data_by_gamma[current_gamma]["S_ei"].append(float(parts[3]))

    return {
        gamma: {
            "q": np.array(block["q"]),
            "S_ee": np.array(block["S_ee"]),
            "S_ii": np.array(block["S_ii"]),
            "S_ei": np.array(block["S_ei"]),
        }
        for gamma, block in data_by_gamma.items()
    }


def render_channel_plot(
    data: dict[float, dict[str, np.ndarray]],
    channel: str,
    title: str,
    filename: str,
) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))

    gammas = sorted(data.keys())
    colors = sns.color_palette("flare", n_colors=len(gammas))

    for idx, gamma in enumerate(gammas):
        q = data[gamma]["q"]
        s_val = data[gamma][channel]
        ax.plot(
            q,
            s_val,
            color=colors[idx],
            linewidth=2.0,
            label=rf"$\Gamma = {gamma:g}$",
        )

    ax.set_title(title, fontweight="bold")
    ax.set_xlim(0, float(np.max(data[gammas[0]]["q"])))

    if channel in ("S_ee", "S_ii"):
        ax.axhline(1.0, color="0.35", linestyle="--", linewidth=1.0)

    ax.set_ylim(bottom=0.0)

    ax.set_xlabel(r"Reduced wavevector $\bar{q} = q/k_F$")
    ax.set_ylabel(r"Static Structure Factor $S(\bar{q})$")
    legend_loc = "upper left" if channel == "S_ii" else "upper right"
    ax.legend(loc=legend_loc, frameon=True, edgecolor="0.7")

    output_path_saved = save_figure(fig, filename)
    plt.close(fig)
    print(f"Saved {output_path_saved}")


def main() -> None:
    configure_seaborn()
    data = parse_gamma_data(DATA_FILE)

    if not data:
        raise SystemExit(
            "Error: No data parsed. Ensure C++ engine is updated and ran --gamma-sweep."
        )

    render_channel_plot(
        data,
        "S_ee",
        r"$S_{ee}(\bar{q})$ — Electron-Electron",
        "Sq_gamma_ee",
    )
    render_channel_plot(
        data,
        "S_ii",
        r"$S_{ii}(\bar{q})$ — Ion-Ion",
        "Sq_gamma_ii",
    )
    render_channel_plot(
        data,
        "S_ei",
        r"$S_{ei}(\bar{q})$ — Electron-Ion Cross",
        "Sq_gamma_ei",
    )


if __name__ == "__main__":
    main()
