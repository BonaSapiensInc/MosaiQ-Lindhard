#!/usr/bin/env python3
"""Render S_ee, S_ii, S_ei(q, omega) contour and 3D surface figures from MosaiQ-Lindhard CLI output."""

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
            colors="black",
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
            colors="black",
            linewidths=0.45,
            alpha=0.75,
        )

    frame_contour_axis(ax)
    ax.set_xlabel(r"Wave vector $q \ [k_F]$")
    ax.set_ylabel(r"Frequency $\omega \ [E_F/\hbar]$")
    ax.set_title(f"{title} — Contour Map", fontweight="bold")
    fig.colorbar(scalar_map, ax=ax, fraction=0.046, pad=0.04, label=label)

    output_path_saved = save_figure(fig, output_path.stem)
    plt.close(fig)
    print(f"Saved {output_path_saved}")


def render_surface(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    plot_grid: np.ma.MaskedArray,
    color_norm: LogNorm | Normalize,
    label: str,
    title: str,
    output_path: Path,
) -> None:
    cmap = cm.rainbow
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
    render_contour(
        q_grid,
        w_grid,
        plot_grid,
        color_norm,
        label,
        title,
        output_path(f"{stem}_contour"),
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

    for stem, col, title in PLOT_SPECS:
        q_grid, w_grid, value_grid = load_gridded(data, col)
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
