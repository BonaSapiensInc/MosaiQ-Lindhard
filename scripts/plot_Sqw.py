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

"""Render S_ee, S_ii, S_ei(q, omega) contour and 3D surface figures for Figure 3."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize

from plot_common import OUTPUT_DIR, PDF_RCPARAMS, output_path, save_figure

DATA_FILE = OUTPUT_DIR / "output_structure_factor.dat"

# Figure 4 background helper bounds (plasmon sector).
FIG4_Q_MIN = 0.1
FIG4_Q_MAX = 4.0
FIG4_OMEGA_YMAX = 10.0

# Electron-channel Figure 3 panels (plasmon / screening sector).
ELECTRON_CHANNEL_Q_MAX = 4.0

# Ion-acoustic contour: macroscopic q extent; omega uses the same simulator grid as electron channels.
S_II_CONTOUR_Q_MAX = 50.0

# Log-spaced isoline styling (white on rainbow fill).
CONTOUR_ISOLINE_COLOR = "white"
CONTOUR_ISOLINE_LINEWIDTH = 1.35
CONTOUR_ISOLINE_ALPHA = 0.88
# Cross-channel panels: sparse signal on a wide log mesh — bold isolines in the active band only.
CROSS_CHANNEL_ISOLINE_LINEWIDTH = 9.0
CROSS_CHANNEL_ISOLINE_ALPHA = 1.0
CROSS_CHANNEL_HALO_COLOR = "#0a0412"
CROSS_CHANNEL_HALO_LINEWIDTH = 22.0
CROSS_CHANNEL_CONTOUR_N_LINES = 12
CROSS_CHANNEL_LINE_VMAX_FLOOR = 0.006
CROSS_CHANNEL_MESH_REFINE = 4

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


def load_gridded(data: np.ndarray, value_col: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_1d = data[:, 0]
    w_1d = data[:, 1]
    values = data[:, value_col]

    q_unique = np.unique(q_1d)
    w_unique = np.unique(w_1d)
    q_grid, w_grid = np.meshgrid(q_unique, w_unique)
    value_grid = values.reshape(len(q_unique), len(w_unique)).T
    return q_grid, w_grid, value_grid


def build_contour_grid(
    data: np.ndarray,
    value_col: int,
    q_min: float,
    q_max: float,
    w_max: float,
    n_omega: int = 400,
    w_min: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate samples onto a rectangular mesh (Figure 4 background fallback)."""
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
    if w_min is not None:
        in_range &= w_col >= w_min - 1.0e-9

    positive_w = w_col[in_range & (w_col > 0.0)]
    if w_min is not None:
        w_axis_lo = w_min
    else:
        w_axis_lo = float(np.min(positive_w)) if positive_w.size else 1.0e-2
        w_axis_lo = max(w_axis_lo, 1.0e-2)
    w_axis = np.geomspace(w_axis_lo, w_max, n_omega)

    value_grid = np.full((n_omega, q_points.size), np.nan)
    for j, q in enumerate(q_points):
        sel = np.isclose(q_col, q, rtol=0.0, atol=1.0e-9) & (w_col <= w_max + 1.0e-9)
        if w_min is not None:
            sel &= w_col >= w_min - 1.0e-9
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


def load_gridded_for_fig4(
    data: np.ndarray, value_col: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_max = min(float(np.max(data[:, 0])), FIG4_Q_MAX)
    w_max = min(float(np.max(data[:, 1])), FIG4_OMEGA_YMAX)
    return build_contour_grid(data, value_col, FIG4_Q_MIN, q_max, w_max, n_omega=250)


def load_channel_grid(
    data: np.ndarray, value_col: int, stem: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the plotting mesh with channel-specific q windows; omega grid is shared."""
    if stem == "S_ii":
        q_cap = min(S_II_CONTOUR_Q_MAX, float(np.max(data[:, 0])))
        channel_rows = data[data[:, 0] <= q_cap + 1.0e-9]
        return load_gridded(channel_rows, value_col)

    q_cap = min(ELECTRON_CHANNEL_Q_MAX, float(np.max(data[:, 0])))
    channel_rows = data[data[:, 0] <= q_cap + 1.0e-9]
    return load_gridded(channel_rows, value_col)


def prepare_plot_grid(grid: np.ndarray) -> tuple[np.ndarray, LogNorm | None]:
    """Map |S| onto a log color scale; clamp vacuum (zero/invalid) to vmin, not white."""
    magnitude = np.abs(grid).astype(float)
    positive = magnitude[np.isfinite(magnitude) & (magnitude > 0.0)]
    norm = choose_log_norm(positive)
    if norm is None:
        return np.nan_to_num(magnitude, nan=0.0), None

    vmin = float(norm.vmin)
    plot_grid = np.where(
        np.isfinite(magnitude) & (magnitude > 0.0),
        np.maximum(magnitude, vmin),
        vmin,
    )
    return plot_grid, norm


def rainbow_log_cmap() -> cm.Colormap:
    """Rainbow with masked/underflow pinned to the purple log-floor."""
    cmap = cm.rainbow.copy()
    floor = cmap(0.0)
    cmap.set_under(floor)
    cmap.set_bad(floor)
    return cmap


def filter_contour_line_levels(
    line_levels: np.ndarray,
    norm: LogNorm,
    *,
    min_vmax_fraction: float = 1.0e-3,
) -> np.ndarray:
    """Drop isolines in the deep-vacuum decade so white lines do not wash out purple fill."""
    vmin = float(norm.vmin)
    vmax = float(norm.vmax)
    floor = max(vmin * 10.0, vmax * min_vmax_fraction)
    return line_levels[line_levels >= floor]


def cross_channel_contour_levels(
    plot_grid: np.ndarray,
    norm: LogNorm,
    *,
    n_lines: int = CROSS_CHANNEL_CONTOUR_N_LINES,
    vmax_floor_fraction: float = CROSS_CHANNEL_LINE_VMAX_FLOOR,
) -> np.ndarray:
    """Place fewer, bolder isolines only where cross-channel spectral weight is non-negligible."""
    vmax = float(norm.vmax)
    lo = max(vmax * vmax_floor_fraction, float(norm.vmin))
    active = plot_grid[np.isfinite(plot_grid) & (plot_grid >= lo)]
    if active.size:
        lo = max(lo, float(np.percentile(active, 8.0)))
    hi = vmax
    if hi <= lo * (1.0 + 1.0e-12):
        return np.array([hi])
    return np.geomspace(lo, hi, n_lines)


def refine_plot_mesh(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    plot_grid: np.ndarray,
    *,
    q_factor: int = CROSS_CHANNEL_MESH_REFINE,
    w_factor: int = CROSS_CHANNEL_MESH_REFINE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Upsample the causal field so cross-channel isolines trace smooth, continuous curves."""
    q_axis = q_grid[0, :]
    w_axis = w_grid[:, 0]
    q_fine = np.linspace(float(q_axis[0]), float(q_axis[-1]), (len(q_axis) - 1) * q_factor + 1)
    w_fine = np.linspace(float(w_axis[0]), float(w_axis[-1]), (len(w_axis) - 1) * w_factor + 1)

    along_q = np.vstack([np.interp(q_fine, q_axis, plot_grid[i, :]) for i in range(len(w_axis))])
    along_w = np.vstack([np.interp(w_fine, w_axis, along_q[:, j]) for j in range(len(q_fine))]).T
    q_fine_grid, w_fine_grid = np.meshgrid(q_fine, w_fine)
    return q_fine_grid, w_fine_grid, along_w


def draw_contour_isolines(
    ax: plt.Axes,
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    plot_grid: np.ndarray,
    color_norm: LogNorm | Normalize,
    line_levels: np.ndarray,
    *,
    cross_channel: bool = False,
    line_width: float | None = None,
    line_alpha: float | None = None,
) -> None:
    width = (
        CROSS_CHANNEL_ISOLINE_LINEWIDTH
        if cross_channel
        else CONTOUR_ISOLINE_LINEWIDTH
        if line_width is None
        else line_width
    )
    alpha = (
        CROSS_CHANNEL_ISOLINE_ALPHA
        if cross_channel
        else CONTOUR_ISOLINE_ALPHA
        if line_alpha is None
        else line_alpha
    )

    contour_q, contour_w, contour_z = q_grid, w_grid, plot_grid
    if cross_channel:
        contour_q, contour_w, contour_z = refine_plot_mesh(q_grid, w_grid, plot_grid)

    contour_set = ax.contour(
        contour_q,
        contour_w,
        contour_z,
        levels=line_levels,
        norm=color_norm,
        colors=CONTOUR_ISOLINE_COLOR,
        linewidths=width,
        alpha=alpha,
        zorder=5,
    )
    if cross_channel:
        contour_set.set(
            path_effects=[
                pe.withStroke(
                    linewidth=CROSS_CHANNEL_HALO_LINEWIDTH,
                    foreground=CROSS_CHANNEL_HALO_COLOR,
                ),
                pe.Normal(),
            ]
        )


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
    ax.set_frame_on(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
    ax.tick_params(top=True, right=True, which="both")


def log_contour_levels(norm: LogNorm, n_filled: int, n_lines: int) -> tuple[np.ndarray, np.ndarray]:
    vmin = float(norm.vmin)
    vmax = float(norm.vmax)
    filled = np.geomspace(vmin, vmax, n_filled)
    lines = np.geomspace(vmin, vmax, n_lines)
    return filled, lines


def render_contour(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    plot_grid: np.ndarray,
    color_norm: LogNorm | Normalize,
    label: str,
    title: str,
    output_path: Path,
    q_xmax: float | None = None,
    omega_ymin: float | None = None,
    omega_ymax: float | None = None,
    contour_linewidth: float | None = None,
    contour_profile: str = "default",
) -> None:
    cross_channel = contour_profile == "cross"
    line_width = (
        CROSS_CHANNEL_ISOLINE_LINEWIDTH
        if cross_channel
        else CONTOUR_ISOLINE_LINEWIDTH
        if contour_linewidth is None
        else contour_linewidth
    )
    line_alpha = CROSS_CHANNEL_ISOLINE_ALPHA if cross_channel else CONTOUR_ISOLINE_ALPHA
    cmap = rainbow_log_cmap()
    scalar_map = ScalarMappable(cmap=cmap, norm=color_norm)
    scalar_map.set_array([])

    fig, ax = plt.subplots(figsize=(10, 7), layout="constrained")
    if isinstance(color_norm, LogNorm):
        fill_levels, line_levels = log_contour_levels(color_norm, n_filled=200, n_lines=72)
        if cross_channel:
            line_levels = cross_channel_contour_levels(plot_grid, color_norm)
        else:
            line_levels = filter_contour_line_levels(line_levels, color_norm)
        ax.contourf(
            q_grid,
            w_grid,
            plot_grid,
            levels=fill_levels,
            cmap=cmap,
            norm=color_norm,
            extend="min",
        )
        draw_contour_isolines(
            ax,
            q_grid,
            w_grid,
            plot_grid,
            color_norm,
            line_levels,
            cross_channel=cross_channel,
            line_width=line_width,
            line_alpha=line_alpha,
        )
    else:
        ax.contourf(q_grid, w_grid, plot_grid, levels=150, cmap=cmap, norm=color_norm)
        draw_contour_isolines(
            ax,
            q_grid,
            w_grid,
            plot_grid,
            color_norm,
            np.geomspace(
                float(color_norm.vmin) if hasattr(color_norm, "vmin") else 0.0,
                float(color_norm.vmax) if hasattr(color_norm, "vmax") else 1.0,
                60,
            ),
            cross_channel=cross_channel,
            line_width=line_width,
            line_alpha=line_alpha,
        )

    frame_contour_axis(ax)
    ax.set_xlabel(r"Wave vector $q \ [k_F]$")
    ax.set_ylabel(r"Frequency $\omega \ [E_F/\hbar]$")
    ax.set_title(f"{title} — Contour Map", fontweight="bold")
    if q_xmax is not None:
        ax.set_xlim(0.0, q_xmax)
    if omega_ymin is not None or omega_ymax is not None:
        y_lo, y_hi = ax.get_ylim()
        ax.set_ylim(
            omega_ymin if omega_ymin is not None else y_lo,
            omega_ymax if omega_ymax is not None else y_hi,
        )
    fig.colorbar(scalar_map, ax=ax, fraction=0.046, pad=0.04, label=label)

    output_path_saved = save_figure(fig, output_path.stem)
    plt.close(fig)
    print(f"Saved {output_path_saved}")


def render_surface(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    plot_grid: np.ndarray,
    color_norm: LogNorm | Normalize,
    label: str,
    title: str,
    output_path: Path,
) -> None:
    cmap = rainbow_log_cmap()
    scalar_map = ScalarMappable(cmap=cmap, norm=color_norm)
    scalar_map.set_array([])

    fig = plt.figure(figsize=(10, 8), layout="constrained")
    ax = fig.add_subplot(111, projection="3d")
    z_surface = plot_grid.astype(float)
    ax.plot_surface(
        q_grid,
        w_grid,
        z_surface,
        cmap=cmap,
        norm=color_norm,
        shade=False,
        edgecolor="#333333",
        linewidth=0.2,
        antialiased=True,
        alpha=1.0,
    )
    ax.set_xlabel(r"$q \ [k_F]$", labelpad=10)
    ax.set_ylabel(r"$\omega \ [E_F/\hbar]$", labelpad=10)
    ax.set_zlabel(label, labelpad=10)
    ax.set_title(f"{title} — 3D Surface", fontweight="bold")
    ax.view_init(elev=35, azim=-120)
    fig.colorbar(scalar_map, ax=ax, fraction=0.04, pad=0.08, label=label)

    output_path_saved = save_figure(fig, output_path.stem)
    plt.close(fig)
    print(f"Saved {output_path_saved}")


def render_channel(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    grid: np.ndarray,
    stem: str,
    label: str,
    title: str,
) -> None:
    plot_grid, norm = prepare_plot_grid(grid)
    color_norm: LogNorm | Normalize = norm if norm is not None else Normalize()

    cross_profile = "cross" if stem == "S_ei" else "default"

    if stem == "S_ii":
        q_xmax = min(S_II_CONTOUR_Q_MAX, float(np.max(q_grid)))
        render_contour(
            q_grid,
            w_grid,
            plot_grid,
            color_norm,
            label,
            title,
            output_path(f"{stem}_contour"),
            q_xmax=q_xmax,
            contour_profile=cross_profile,
        )
        # 3D companion stays on the low-q electron-sector mesh for readability.
        electron_rows_mask = q_grid[0, :] <= ELECTRON_CHANNEL_Q_MAX + 1.0e-9
        surface_q = q_grid[:, electron_rows_mask]
        surface_w = w_grid[:, electron_rows_mask]
        surface_grid = grid[:, electron_rows_mask]
        surface_plot, surface_norm = prepare_plot_grid(surface_grid)
        surface_color: LogNorm | Normalize = (
            surface_norm if surface_norm is not None else Normalize()
        )
        render_surface(
            surface_q,
            surface_w,
            surface_plot,
            surface_color,
            label,
            title,
            output_path(f"{stem}_3d"),
        )
        return

    render_contour(
        q_grid,
        w_grid,
        plot_grid,
        color_norm,
        label,
        title,
        output_path(f"{stem}_contour"),
        contour_profile=cross_profile,
    )
    render_surface(
        q_grid,
        w_grid,
        plot_grid,
        color_norm,
        label,
        title,
        output_path(f"{stem}_3d"),
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

    q_data_max = float(np.max(data[:, 0]))
    if q_data_max + 1.0e-9 < S_II_CONTOUR_Q_MAX:
        print(
            f"Warning: data q_max={q_data_max:g} < {S_II_CONTOUR_Q_MAX:g}; "
            "re-run mosaiq_simulator to refresh the wide S_ii mesh."
        )

    for stem, col, title in PLOT_SPECS:
        q_grid, w_grid, value_grid = load_channel_grid(data, col, stem)
        render_channel(
            q_grid,
            w_grid,
            value_grid,
            stem,
            label=title.split("—")[0].strip(),
            title=title,
        )


if __name__ == "__main__":
    main()
