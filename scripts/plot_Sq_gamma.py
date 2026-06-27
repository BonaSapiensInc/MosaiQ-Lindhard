#!/usr/bin/env python3

# ==========================================================================
# MOSAIQ-LINDHARD ENGINE
# Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
# * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
# Free for non-commercial personal and academic research use only.
# Commercial, governmental, and public institutional use requires a
# separate paid license. See LICENSE file for details.
# * Contact: kim.ingee@bonasapiens.com
# ==========================================================================

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
INSET_Q_MAX = 4.0
S_EE_INSET_Q_MIN = 1.4
S_EE_INSET_Q_MAX = 2.0
# Axes-fraction box spanning main-axis q ≈ 10–35 (xlim 0–50 → 0.20–0.70).
INSET_BBOX = (0.20, 0.48, 0.50, 0.38)
# S_ee inset: anchor top-left, expand down/right while clearing lower-right legend.
S_EE_INSET_TOP_LEFT = (0.20, 0.86)
S_EE_INSET_SIZE = (0.56, 0.44)
S_EI_YSCALE = 1.0e3


def channel_yvalues(channel: str, values: np.ndarray) -> np.ndarray:
    if channel == "S_ei":
        return values * S_EI_YSCALE
    return values


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


def inset_bbox_from_topleft(
    left: float, top: float, width: float, height: float
) -> tuple[float, float, float, float]:
    return (left, top - height, width, height)


def add_low_q_inset(
    ax: plt.Axes,
    data: dict[float, dict[str, np.ndarray]],
    gammas: list[float],
    colors: list,
    channel: str,
    *,
    left: float | None = None,
    bottom: float | None = None,
    width: float | None = None,
    height: float | None = None,
) -> None:
    default_left, default_bottom, default_width, default_height = INSET_BBOX
    axins = ax.inset_axes(
        [
            left if left is not None else default_left,
            bottom if bottom is not None else default_bottom,
            width if width is not None else default_width,
            height if height is not None else default_height,
        ]
    )
    axins.set_facecolor("white")

    inset_q_min = S_EE_INSET_Q_MIN if channel == "S_ee" else 0.0
    inset_q_max = S_EE_INSET_Q_MAX if channel == "S_ee" else INSET_Q_MAX
    inset_vals: list[float] = []
    for idx, gamma in enumerate(gammas):
        q = data[gamma]["q"]
        s_val = channel_yvalues(channel, data[gamma][channel])
        mask = (q >= inset_q_min) & (q <= inset_q_max)
        axins.plot(
            q[mask],
            s_val[mask],
            color=colors[idx],
            linewidth=1.5,
        )
        inset_vals.extend(s_val[mask].tolist())

    axins.set_xlim(inset_q_min, inset_q_max)
    if channel == "S_ee":
        axins.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7, zorder=0)
        axins.set_ylim(0.0, 0.5)
    elif channel == "S_ii":
        axins.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7, zorder=0)
        axins.set_ylim(0.0, max(inset_vals) * 1.08)
    else:
        axins.axhline(0.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7, zorder=0)
        inset_max = max(inset_vals)
        axins.set_ylim(-0.02 * inset_max, inset_max * 1.15)

    axins.tick_params(labelsize=7)
    axins.grid(True, alpha=0.25)
    for spine in axins.spines.values():
        spine.set_linewidth(0.8)


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
        s_val = channel_yvalues(channel, data[gamma][channel])
        ax.plot(
            q,
            s_val,
            color=colors[idx],
            linewidth=2.0,
            label=rf"$\Gamma = {gamma:g}$",
        )

    ax.set_title(title, fontweight="bold")
    ax.set_xlim(0.0, 50.0)

    if channel in ("S_ee", "S_ii"):
        ax.axhline(
            1.0,
            color="black",
            linestyle="--",
            linewidth=1.0,
            alpha=0.7,
            zorder=0,
        )
        channel_max = max(float(np.max(data[g][channel])) for g in gammas)
        ax.set_ylim(0.0, max(1.08, channel_max * 1.05))
    elif channel == "S_ei":
        ax.axhline(
            0.0,
            color="black",
            linestyle="--",
            linewidth=1.0,
            alpha=0.7,
            zorder=0,
        )
        channel_max = max(float(np.max(channel_yvalues(channel, data[g][channel]))) for g in gammas)
        ax.set_ylim(-0.02 * channel_max, channel_max * 1.15)

    ax.set_xlabel(r"Reduced wavevector $\bar{q} = q/k_F$")
    if channel == "S_ei":
        ax.set_ylabel(r"Static Structure Factor $S(\bar{q})$ ($\times 10^{-3}$)")
    else:
        ax.set_ylabel(r"Static Structure Factor $S(\bar{q})$")
    legend_loc = "lower right" if channel in ("S_ee", "S_ii") else "upper right"
    ax.legend(loc=legend_loc, frameon=True, edgecolor="0.7")

    if channel == "S_ii":
        inset_left, inset_bottom = 0.42, 0.38
        inset_width, inset_height = None, None
    elif channel == "S_ee":
        ee_left, ee_top = S_EE_INSET_TOP_LEFT
        ee_width, ee_height = S_EE_INSET_SIZE
        inset_left, inset_bottom, inset_width, inset_height = inset_bbox_from_topleft(
            ee_left, ee_top, ee_width, ee_height
        )
    else:
        inset_left, inset_bottom = 0.42, 0.22
        inset_width, inset_height = None, None

    add_low_q_inset(
        ax,
        data,
        gammas,
        colors,
        channel,
        left=inset_left,
        bottom=inset_bottom,
        width=inset_width,
        height=inset_height,
    )

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
