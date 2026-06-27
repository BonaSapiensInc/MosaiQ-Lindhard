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

"""Render S_ee, S_ii, S_ei(q, omega) focused and wide-range contour figures from MosaiQ-Lindhard CLI output."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize

from plot_common import OUTPUT_DIR, PDF_RCPARAMS, output_path, save_figure

DATA_FILE = OUTPUT_DIR / "output_structure_factor.dat"

MICRO_Q_MIN = 0.1
MICRO_Q_MAX = 4.0
MACRO_Q_MIN = 0.1
MACRO_Q_MAX = 50.0

FOCUSED_OMEGA_YMAX = 10.0
MACRO_OMEGA_YMAX = MACRO_Q_MAX * MACRO_Q_MAX

PLOT_SPECS: tuple[tuple[str, int, str], ...] = (
    ("S_ee", 3, r"$S_{ee}(q, \omega)$ — electron-electron"),
    ("S_ii", 5, r"$S_{ii}(q, \omega)$ — ion-ion (ion-acoustic)"),
    ("S_ei", 7, r"$S_{ei}(q, \omega)$ — electron-ion cross"),
)


def configure_seaborn() -> None:
    sns.set_theme(
        style="ticks",
        context="paper",
        font_scale=1.05,
        rc={
            "font.family": "serif",
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "figure.dpi": 100,
            **PDF_RCPARAMS,
        },
    )


def build_contour_grid(
    data: np.ndarray,
    value_col: int,
    q_min: float,
    q_max: float,
    w_max: float,
    n_omega: int = 400,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate adaptive (q, omega) samples onto a rectangular mesh for contour plotting."""
    q_col = data[:, 0]
    w_col = data[:, 1]
    v_col = data[:, value_col]

    q_points = np.unique(q_col[(q_col >= q_min - 1.0e-9) & (q_col <= q_max + 1.0e-9)])
    if q_points.size == 0:
        raise ValueError(f"No q samples in [{q_min}, {q_max}]")

    in_range = (
        (q_col >= q_min - 1.0e-9)
        & (q_col <= q_max + 1.0e-9)
        & (w_col <= w_max + 1.0e-9)
    )
    positive_w = w_col[in_range & (w_col > 0.0)]
    w_min = float(np.min(positive_w)) if positive_w.size else 1.0e-2
    w_min = max(w_min, 1.0e-2)
    w_axis = np.geomspace(w_min, w_max, n_omega)

    value_grid = np.full((n_omega, q_points.size), np.nan)
    for j, q in enumerate(q_points):
        sel = np.isclose(q_col, q, rtol=0.0, atol=1.0e-9) & (w_col <= w_max + 1.0e-9)
        if not np.any(sel):
            continue
        w_slice = w_col[sel]
        v_slice = v_col[sel]
        order = np.argsort(w_slice)
        w_slice = w_slice[order]
        v_slice = v_slice[order]
        value_grid[:, j] = np.interp(
            w_axis, w_slice, v_slice, left=np.nan, right=np.nan
        )

    q_grid, w_grid = np.meshgrid(q_points, w_axis)
    return q_grid, w_grid, value_grid


def load_gridded(data: np.ndarray, value_col: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a contour mesh over the full adaptive dataset (used by Figure 4 background)."""
    q_max = float(np.max(data[:, 0]))
    w_max = float(np.max(data[:, 1]))
    return build_contour_grid(data, value_col, MICRO_Q_MIN, q_max, w_max, n_omega=400)


def subset_gridded(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    value_grid: np.ndarray,
    q_min: float,
    q_max: float,
    w_max: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_axis = q_grid[0, :]
    w_axis = w_grid[:, 0]
    q_sel = (q_axis >= q_min - 1.0e-9) & (q_axis <= q_max + 1.0e-9)
    w_sel = np.ones_like(w_axis, dtype=bool)
    if w_max is not None:
        w_sel = w_axis <= w_max + 1.0e-9

    idx = np.ix_(w_sel, q_sel)
    return q_grid[idx], w_grid[idx], value_grid[idx]


def prepare_plot_grid(grid: np.ndarray) -> tuple[np.ma.MaskedArray, LogNorm | None]:
    plot_grid = np.ma.masked_invalid(np.abs(grid))
    norm = choose_log_norm(plot_grid.filled(np.nan))
    if norm is not None:
        plot_grid = np.ma.masked_less_equal(plot_grid, 0.0)
    return plot_grid, norm


def choose_log_norm(grid: np.ndarray) -> LogNorm | None:
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        return None

    magnitude = np.abs(finite)
    abs_max = float(np.max(magnitude))
    if abs_max == 0.0:
        return None

    positive = magnitude[magnitude > 0.0]
    min_pos = float(np.min(positive)) if positive.size else abs_max * 1.0e-6
    vmin = max(min_pos, abs_max * 1.0e-6, 1.0e-12)
    return LogNorm(vmin=vmin, vmax=abs_max)


def frame_contour_axis(ax: plt.Axes) -> None:
    """Closed rectangular frame for the 2D contour panel."""
    ax.set_frame_on(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
    ax.tick_params(top=True, right=True, which="both")


def log_contour_levels(norm: LogNorm, n_filled: int, n_lines: int) -> tuple[np.ndarray, np.ndarray]:
    """Log-spaced levels so contour density is uniform across decades of magnitude."""
    vmin = float(norm.vmin)
    vmax = float(norm.vmax)
    filled = np.geomspace(vmin, vmax, n_filled)
    lines = np.geomspace(vmin, vmax, n_lines)
    return filled, lines


def render_contour(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    plot_grid: np.ma.MaskedArray,
    color_norm: LogNorm | Normalize,
    label: str,
    title: str,
    output_path: Path,
    omega_ymax: float,
) -> None:
    cmap = cm.rainbow
    scalar_map = ScalarMappable(cmap=cmap, norm=color_norm)
    scalar_map.set_array([])

    fig, ax = plt.subplots(figsize=(10, 7), layout="constrained")
    if isinstance(color_norm, LogNorm):
        fill_levels, line_levels = log_contour_levels(color_norm, n_filled=200, n_lines=72)
        ax.contourf(q_grid, w_grid, plot_grid, levels=fill_levels, cmap=cmap, norm=color_norm)
        ax.contour(
            q_grid,
            w_grid,
            plot_grid,
            levels=line_levels,
            norm=color_norm,
            colors="white",
            linewidths=0.45,
            alpha=0.75,
        )
    else:
        ax.contourf(q_grid, w_grid, plot_grid, levels=150, cmap=cmap, norm=color_norm)
        ax.contour(
            q_grid,
            w_grid,
            plot_grid,
            levels=60,
            norm=color_norm,
            colors="white",
            linewidths=0.45,
            alpha=0.75,
        )

    frame_contour_axis(ax)
    ax.set_xlabel(r"Reduced wavevector $\bar{q} = q/k_\mathrm{F}$")
    ax.set_ylabel(r"Reduced frequency $\bar{\omega} = \hbar\omega/\epsilon_\mathrm{F}$")
    ax.set_ylim(0.0, omega_ymax)
    ax.set_title(title, fontweight="bold")
    fig.colorbar(scalar_map, ax=ax, fraction=0.046, pad=0.04, label=label)

    output_path_saved = save_figure(fig, output_path.stem)
    plt.close(fig)
    print(f"Saved {output_path_saved}")


def render_channel(
    data: np.ndarray,
    value_col: int,
    stem: str,
    label: str,
    title: str,
) -> None:
    micro_q, micro_w, micro_grid = build_contour_grid(
        data, value_col, MICRO_Q_MIN, MICRO_Q_MAX, w_max=FOCUSED_OMEGA_YMAX, n_omega=200
    )
    macro_w_max = MACRO_OMEGA_YMAX
    macro_q, macro_w, macro_grid = build_contour_grid(
        data, value_col, MACRO_Q_MIN, MACRO_Q_MAX, w_max=macro_w_max, n_omega=400
    )

    micro_plot, micro_norm = prepare_plot_grid(micro_grid)
    macro_plot, macro_norm = prepare_plot_grid(macro_grid)

    micro_color_norm: LogNorm | Normalize = (
        micro_norm if micro_norm is not None else Normalize()
    )
    macro_color_norm: LogNorm | Normalize = (
        macro_norm if macro_norm is not None else Normalize()
    )

    render_contour(
        micro_q,
        micro_w,
        micro_plot,
        micro_color_norm,
        label,
        f"{title} — Focused Contour ($\\bar{{q}} \\leq {MICRO_Q_MAX:.1f}$)",
        output_path(f"{stem}_contour"),
        FOCUSED_OMEGA_YMAX,
    )
    render_contour(
        macro_q,
        macro_w,
        macro_plot,
        macro_color_norm,
        label,
        f"{title} — Extended Contour ($\\bar{{q}} \\leq {MACRO_Q_MAX:.1f}$)",
        output_path(f"{stem}_wide_contour"),
        MACRO_OMEGA_YMAX,
    )


def main() -> None:
    configure_seaborn()

    if not DATA_FILE.is_file():
        raise SystemExit(f"Error: data file not found: {DATA_FILE}")

    data = np.loadtxt(DATA_FILE, comments="#")
    if data.ndim != 2 or data.shape[1] < 8:
        raise SystemExit(
            f"Error: expected at least 8 columns in {DATA_FILE}, got shape {data.shape}"
        )

    for stem, col, title in PLOT_SPECS:
        render_channel(
            data,
            col,
            stem,
            label=title.split("—")[0].strip(),
            title=title,
        )


if __name__ == "__main__":
    main()
