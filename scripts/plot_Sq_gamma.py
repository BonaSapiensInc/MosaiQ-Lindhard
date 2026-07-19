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

"""Render Static Structure Factors S(q) per channel with Gamma sweeps (PRE style).

Phase Z4: frames multi-component Zeta-RPA static S(q) with emphasis on
extreme-coupling (Γ ≥ 100) correlation-peak stabilization.
"""

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
# Ion inset must reach the extreme-Γ coordination spike (~q̄ ≈ 5–6).
S_II_INSET_Q_MAX = 8.0
# Default crystallization window; refined from data when peaks are locatable.
S_EE_INSET_Q_MIN = 1.15
S_EE_INSET_Q_MAX = 2.20
# Axes-fraction box spanning main-axis q ≈ 10–35 (xlim 0–50 → 0.20–0.70).
INSET_BBOX = (0.20, 0.48, 0.50, 0.38)
# S_ee inset: anchor top-left, expand down/right while clearing lower-right legend.
S_EE_INSET_TOP_LEFT = (0.20, 0.86)
S_EE_INSET_SIZE = (0.56, 0.44)
# Cross channel is O(10^{-3})–O(10^{-6}); keep 10^3 scale for Γ=50 drag peak.
S_EI_YSCALE = 1.0e3
# Γ ≥ this value is treated as Zeta-dominated extreme coupling.
EXTREME_GAMMA_FLOOR = 100.0

# Highest Gamma gets the plainest linestyle; weaker coupling adds dash complexity.
GAMMA_LINE_STYLES: tuple[str | tuple[int, tuple[int, ...]], ...] = (
    "solid",
    "dashed",
    (0, (4, 2, 1, 2)),  # dot-dashed
    (0, (1, 2, 3, 2, 1, 2)),  # dot-dot-dashed
)
# Rank 0 = highest Γ: thickest stroke for stabilized Zeta structure.
GAMMA_LINE_WIDTHS: tuple[float, ...] = (2.8, 2.3, 1.7, 1.35)


def gamma_linestyle(gamma: float, gammas: list[float]) -> str | tuple[int, tuple[int, ...]]:
    rank = sorted(gammas, reverse=True).index(gamma)
    return GAMMA_LINE_STYLES[min(rank, len(GAMMA_LINE_STYLES) - 1)]


def gamma_linewidth(gamma: float, gammas: list[float]) -> float:
    rank = sorted(gammas, reverse=True).index(gamma)
    return GAMMA_LINE_WIDTHS[min(rank, len(GAMMA_LINE_WIDTHS) - 1)]


def gamma_alpha(gamma: float) -> float:
    return 1.0 if gamma >= EXTREME_GAMMA_FLOOR else 0.82


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


def parse_gamma_data(path: Path) -> tuple[dict[float, dict[str, np.ndarray]], str]:
    if not path.is_file():
        raise SystemExit(f"Error: data file not found: {path}")

    data_by_gamma: dict[float, dict[str, list[float]]] = defaultdict(
        lambda: {"q": [], "S_ee": [], "S_ii": [], "S_ei": []}
    )
    current_gamma: float | None = None
    pathway = "ZetaRPA"

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            pathway_match = re.search(r"ResponsePathway\s*=\s*(\S+)", line)
            if pathway_match:
                pathway = pathway_match.group(1)
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

    data = {
        gamma: {
            "q": np.array(block["q"]),
            "S_ee": np.array(block["S_ee"]),
            "S_ii": np.array(block["S_ii"]),
            "S_ei": np.array(block["S_ei"]),
        }
        for gamma, block in data_by_gamma.items()
    }
    return data, pathway


def inset_bbox_from_topleft(
    left: float, top: float, width: float, height: float
) -> tuple[float, float, float, float]:
    return (left, top - height, width, height)


def resolve_ee_inset_window(
    data: dict[float, dict[str, np.ndarray]], gammas: list[float]
) -> tuple[float, float]:
    """Frame the first S_ee coordination sector around q̄ ∈ [1.15, 2.20].

    Extreme-coupling Zeta-RPA often flattens the classical crystallization peak into
    a rising shoulder. Prefer an interior local max when Γ ≥ 100 still shows one;
    otherwise keep the default window so weak- and strong-Γ curves stay co-framed.
    """
    extreme = [g for g in gammas if g >= EXTREME_GAMMA_FLOOR]
    peak_qs: list[float] = []
    for gamma in extreme:
        q = data[gamma]["q"]
        s = data[gamma]["S_ee"]
        mask = (q >= 1.0) & (q <= 2.6)
        if not np.any(mask):
            continue
        qq = q[mask]
        ss = s[mask]
        if len(ss) < 3:
            continue
        for i in range(1, len(ss) - 1):
            if ss[i] > ss[i - 1] and ss[i] >= ss[i + 1] and ss[i] > 0.02:
                peak_qs.append(float(qq[i]))
                break

    if not peak_qs:
        return S_EE_INSET_Q_MIN, S_EE_INSET_Q_MAX

    q_c = float(np.median(peak_qs))
    half = 0.50
    return max(1.0, q_c - half), min(2.6, q_c + half)


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
    ee_q_min: float = S_EE_INSET_Q_MIN,
    ee_q_max: float = S_EE_INSET_Q_MAX,
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

    inset_q_min = ee_q_min if channel == "S_ee" else 0.0
    if channel == "S_ee":
        inset_q_max = ee_q_max
    elif channel == "S_ii":
        inset_q_max = S_II_INSET_Q_MAX
    else:
        inset_q_max = INSET_Q_MAX
    inset_vals: list[float] = []
    for idx, gamma in enumerate(gammas):
        q = data[gamma]["q"]
        s_val = channel_yvalues(channel, data[gamma][channel])
        mask = (q >= inset_q_min) & (q <= inset_q_max)
        axins.plot(
            q[mask],
            s_val[mask],
            color=colors[idx],
            linewidth=gamma_linewidth(gamma, gammas) * 0.72,
            linestyle=gamma_linestyle(gamma, gammas),
            alpha=gamma_alpha(gamma),
            solid_capstyle="round",
        )
        inset_vals.extend(s_val[mask].tolist())

    axins.set_xlim(inset_q_min, inset_q_max)
    if not inset_vals:
        return

    inset_min = float(min(inset_vals))
    inset_max = float(max(inset_vals))
    span = max(inset_max - inset_min, 1.0e-12)

    if channel == "S_ee":
        axins.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7, zorder=0)
        # Tight frame around the crystallization peak (no empty 0–0.5 void).
        y0 = max(0.0, inset_min - 0.08 * span)
        y1 = inset_max + 0.18 * span
        axins.set_ylim(y0, y1)
    elif channel == "S_ii":
        axins.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7, zorder=0)
        axins.set_ylim(0.0, inset_max * 1.12)
    else:
        axins.axhline(0.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7, zorder=0)
        axins.set_ylim(inset_min - 0.08 * span, inset_max + 0.18 * span)

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
    ee_q_min, ee_q_max = resolve_ee_inset_window(data, gammas)

    for idx, gamma in enumerate(gammas):
        q = data[gamma]["q"]
        s_val = channel_yvalues(channel, data[gamma][channel])
        ax.plot(
            q,
            s_val,
            color=colors[idx],
            linewidth=gamma_linewidth(gamma, gammas),
            linestyle=gamma_linestyle(gamma, gammas),
            alpha=gamma_alpha(gamma),
            solid_capstyle="round",
            label=rf"$\Gamma = {gamma:g}$",
            zorder=3 + sorted(gammas).index(gamma),
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
        scaled = [channel_yvalues(channel, data[g][channel]) for g in gammas]
        channel_min = min(float(np.min(s)) for s in scaled)
        channel_max = max(float(np.max(s)) for s in scaled)
        span = max(channel_max - channel_min, 1.0e-12)
        ax.set_ylim(channel_min - 0.05 * span, channel_max + 0.12 * span)

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
        ee_q_min=ee_q_min,
        ee_q_max=ee_q_max,
    )

    output_path_saved = save_figure(fig, filename)
    plt.close(fig)
    print(f"Saved {output_path_saved}")
    if channel == "S_ee":
        print(f"  S_ee inset window: q ∈ [{ee_q_min:.3f}, {ee_q_max:.3f}]")


def main() -> None:
    configure_seaborn()
    data, pathway = parse_gamma_data(DATA_FILE)

    if not data:
        raise SystemExit(
            "Error: No data parsed. Ensure C++ engine is updated and ran --gamma-sweep."
        )

    gammas = sorted(data.keys())
    print(f"Pathway={pathway}; Gamma blocks={gammas}")
    if len(gammas) < 4:
        print(
            "Warning: expected Γ ∈ {10,50,100,150}; plot will use whatever blocks are present.",
            file=sys.stderr,
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
