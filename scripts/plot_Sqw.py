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

"""Render S_ee, S_ii, S_ei(q, omega) contour figures for Figure 4."""

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

# Ion-acoustic contour: macroscopic q extent with Golden Scale omega zoom.
S_II_CONTOUR_Q_MAX = 50.0
ION_CHANNEL_OMEGA_MAX = 0.20
ION_CHANNEL_OMEGA_DATA_PAD = 0.02

# Tidy display frames (match t0_analytic contour styling).
DISPLAY_MESH_POINTS = 400
ELECTRON_DISPLAY_BOUNDS = (0.0, 4.0, 0.0, 3.0)
ION_DISPLAY_BOUNDS = (0.0, 50.0, 0.0, ION_CHANNEL_OMEGA_MAX)

# Backward-compatible aliases.
ION_CHANNEL_OMEGA_DATA_MIN = 0.01
ION_CHANNEL_OMEGA_MIN = ION_CHANNEL_OMEGA_DATA_MIN

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
    n_q: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate samples onto a rectangular mesh (Figure 4 background fallback)."""
    q_col = data[:, 0]
    w_col = data[:, 1]
    v_col = data[:, value_col]

    in_q = (q_col >= q_min - 1.0e-9) & (q_col <= q_max + 1.0e-9)
    if n_q is not None:
        q_points = np.linspace(q_min, q_max, n_q)
        source_q = np.unique(q_col[in_q])
    else:
        q_points = np.unique(q_col[in_q])
        source_q = q_points
    if q_points.size == 0 or source_q.size == 0:
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
        q_lookup = q
        if n_q is not None:
            q_lookup = float(source_q[np.argmin(np.abs(source_q - q))])
        sel = np.isclose(q_col, q_lookup, rtol=0.0, atol=1.0e-9) & (w_col <= w_max + 1.0e-9)
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


def crop_omega_band(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    grid: np.ndarray,
    omega_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Restrict contour data to the ion-acoustic frequency window."""
    band = w_grid[:, 0] <= omega_max + 1.0e-9
    return q_grid[band, :], w_grid[band, :], grid[band, :]


def ion_channel_data_omega_cap() -> float:
    """Include one simulator step beyond the display ceiling for edge interpolation."""
    return ION_CHANNEL_OMEGA_MAX + ION_CHANNEL_OMEGA_DATA_PAD


def interpolate_to_display_mesh(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    grid: np.ndarray,
    q_lo: float,
    q_hi: float,
    w_lo: float,
    w_hi: float,
    *,
    n_points: int = DISPLAY_MESH_POINTS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resample onto a uniform rectangle so contours fill the axis frame exactly."""
    q_axis = np.linspace(q_lo, q_hi, n_points)
    w_axis = np.linspace(w_lo, w_hi, n_points)
    q_src = q_grid[0, :]
    w_src = w_grid[:, 0]
    source = grid.astype(float)

    along_q = np.vstack(
        [
            np.interp(q_axis, q_src, source[i, :], left=np.nan, right=np.nan)
            for i in range(len(w_src))
        ]
    )
    along_w = np.vstack(
        [
            np.interp(w_axis, w_src, along_q[:, j], left=np.nan, right=np.nan)
            for j in range(len(q_axis))
        ]
    ).T
    q_display, w_display = np.meshgrid(q_axis, w_axis)
    return q_display, w_display, along_w


def fill_display_boundary(grid: np.ndarray) -> np.ndarray:
    """Extend nearest valid samples across every row/column gap (plasmon backgrounds)."""
    out = np.array(grid, dtype=float, copy=True)
    for i in range(out.shape[0]):
        row = out[i]
        finite = np.isfinite(row)
        if not np.any(finite):
            continue
        j_first = int(np.argmax(finite))
        j_last = len(row) - 1 - int(np.argmax(finite[::-1]))
        row[:j_first] = row[j_first]
        row[j_last + 1 :] = row[j_last]

    for j in range(out.shape[1]):
        col = out[:, j]
        finite = np.isfinite(col)
        if not np.any(finite):
            continue
        i_first = int(np.argmax(finite))
        i_last = len(col) - 1 - int(np.argmax(finite[::-1]))
        col[:i_first] = col[i_first]
        col[i_last + 1 :] = col[i_last]
    return out


def fill_display_border(grid: np.ndarray) -> np.ndarray:
    """Extend nearest interior samples into the outermost display rows/columns only."""
    out = np.array(grid, dtype=float, copy=True)
    nrows, ncols = out.shape
    if nrows < 3 or ncols < 3:
        return out

    def nearest_valid_row(start: int, step: int) -> int | None:
        index = start
        while 0 <= index < nrows:
            if np.any(np.isfinite(out[index])):
                return index
            index += step
        return None

    def nearest_valid_col(start: int, step: int) -> int | None:
        index = start
        while 0 <= index < ncols:
            if np.any(np.isfinite(out[:, index])):
                return index
            index += step
        return None

    top = nearest_valid_row(1, 1)
    if top is not None:
        missing = ~np.isfinite(out[0])
        out[0, missing] = out[top, missing]

    bottom = nearest_valid_row(nrows - 2, -1)
    if bottom is not None:
        missing = ~np.isfinite(out[-1])
        out[-1, missing] = out[bottom, missing]

    left = nearest_valid_col(1, 1)
    if left is not None:
        missing = ~np.isfinite(out[:, 0])
        out[missing, 0] = out[missing, left]

    right = nearest_valid_col(ncols - 2, -1)
    if right is not None:
        missing = ~np.isfinite(out[:, -1])
        out[missing, -1] = out[missing, right]

    return out


def format_contour_display_axes(
    ax: plt.Axes,
    q_lo: float,
    q_hi: float,
    w_lo: float,
    w_hi: float,
) -> None:
    """Pin axis limits and major ticks to clean values like t0_analytic panels."""
    ax.set_xlim(q_lo, q_hi)
    ax.set_ylim(w_lo, w_hi)
    ax.margins(x=0, y=0)
    if q_hi <= 4.0 + 1.0e-9:
        ax.set_xticks(np.arange(int(q_lo), int(q_hi) + 1, 1))
    else:
        ax.set_xticks(np.arange(int(q_lo), int(q_hi) + 1, 10))
    if w_hi <= 0.25 + 1.0e-9:
        ax.set_yticks(np.linspace(w_lo, w_hi, 5))
    else:
        ax.set_yticks(np.arange(int(w_lo), int(w_hi) + 1, 1))


def prepare_ion_channel_mesh(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map ion-acoustic data onto the Golden Scale display rectangle."""
    cropped_q, cropped_w, cropped_grid = crop_omega_band(
        q_grid, w_grid, grid, ion_channel_data_omega_cap()
    )
    return interpolate_to_display_mesh(cropped_q, cropped_w, cropped_grid, *ION_DISPLAY_BOUNDS)


def prepare_ion_plot_grid(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    grid: np.ndarray,
) -> tuple[np.ndarray, LogNorm | None]:
    """Color-normalize on the Golden Scale band only, excluding the padding row."""
    magnitude = np.abs(grid).astype(float)
    in_band = w_grid <= ION_CHANNEL_OMEGA_MAX + 1.0e-9
    band_magnitude = np.where(in_band, magnitude, np.nan)
    positive = band_magnitude[np.isfinite(band_magnitude) & (band_magnitude > 0.0)]
    norm = choose_log_norm(positive)
    if norm is None:
        return np.nan_to_num(magnitude, nan=0.0), None

    vmin = float(norm.vmin)
    plot_grid = np.where(
        in_band & np.isfinite(magnitude) & (magnitude > 0.0),
        np.maximum(magnitude, vmin),
        np.nan,
    )
    return plot_grid, norm


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
    q_bounds: tuple[float, float] | None = None,
    w_bounds: tuple[float, float] | None = None,
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
    if q_bounds is not None and w_bounds is not None:
        format_contour_display_axes(ax, *q_bounds, *w_bounds)
    fig.colorbar(scalar_map, ax=ax, fraction=0.046, pad=0.04, label=label)

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
    cross_profile = "cross" if stem == "S_ei" else "default"
    magnitude = np.abs(grid)

    if stem == "S_ii":
        q_display, w_display, value_display = prepare_ion_channel_mesh(
            q_grid, w_grid, magnitude
        )
        plot_grid, color_norm = prepare_ion_plot_grid(q_display, w_display, value_display)
        if color_norm is None:
            raise RuntimeError("No finite positive amplitudes for S_ii contour")
        q_bounds = ION_DISPLAY_BOUNDS[:2]
        w_bounds = ION_DISPLAY_BOUNDS[2:]
    else:
        q_display, w_display, value_display = interpolate_to_display_mesh(
            q_grid, w_grid, magnitude, *ELECTRON_DISPLAY_BOUNDS
        )
        plot_grid, norm = prepare_plot_grid(value_display)
        color_norm: LogNorm | Normalize = norm if norm is not None else Normalize()
        q_bounds = ELECTRON_DISPLAY_BOUNDS[:2]
        w_bounds = ELECTRON_DISPLAY_BOUNDS[2:]

    render_contour(
        q_display,
        w_display,
        plot_grid,
        color_norm,
        label,
        title,
        output_path(f"{stem}_contour"),
        q_bounds=q_bounds,
        w_bounds=w_bounds,
        contour_profile=cross_profile,
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
