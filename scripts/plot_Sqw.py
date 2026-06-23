#!/usr/bin/env python3
"""Render S_ee, S_ii, S_ei(q, omega) from SOPHUS multi-component RPA output."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILE = PROJECT_ROOT / "output" / "output_structure_factor.dat"
OUTPUT_DIR = PROJECT_ROOT / "output"

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
            "savefig.dpi": 300,
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


def render_component(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    grid: np.ndarray,
    label: str,
    title: str,
    output_path: Path,
) -> None:
    plot_grid = np.ma.masked_invalid(np.abs(grid))
    norm = choose_log_norm(plot_grid.filled(np.nan))
    if norm is not None:
        plot_grid = np.ma.masked_less_equal(plot_grid, 0.0)

    cmap = cm.rainbow
    color_norm: LogNorm | Normalize = norm if norm is not None else Normalize()
    scalar_map = ScalarMappable(cmap=cmap, norm=color_norm)
    scalar_map.set_array([])

    fig = plt.figure(figsize=(12, 10), layout="constrained")

    ax1 = fig.add_subplot(211)
    ax1.contourf(
        q_grid,
        w_grid,
        plot_grid,
        levels=150,
        cmap=cmap,
        norm=color_norm,
    )
    if norm is not None:
        ax1.contour(
            q_grid,
            w_grid,
            plot_grid,
            levels=25,
            norm=color_norm,
            colors="black",
            linewidths=0.8,
            alpha=0.8,
        )
    sns.despine(ax=ax1)
    ax1.set_xlabel(r"Wave vector $q \ [k_F]$")
    ax1.set_ylabel(r"Frequency $\omega \ [E_F/\hbar]$")
    ax1.set_title(f"{title} — Contour Map", fontweight="bold")

    ax2 = fig.add_subplot(212, projection="3d")
    z_surface = plot_grid.astype(float)
    ax2.plot_surface(
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
    ax2.set_xlabel(r"$q \ [k_F]$", labelpad=10)
    ax2.set_ylabel(r"$\omega \ [E_F/\hbar]$", labelpad=10)
    ax2.set_zlabel(label, labelpad=10)
    ax2.set_title(f"{title} — 3D Surface", fontweight="bold")
    ax2.view_init(elev=35, azim=-120)

    fig.colorbar(scalar_map, ax=[ax1, ax2], fraction=0.025, pad=0.04, label=label)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}")


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
        render_component(
            q_grid,
            w_grid,
            value_grid,
            label=title.split("—")[0].strip(),
            title=title,
            output_path=OUTPUT_DIR / f"{stem}_plot.png",
        )


if __name__ == "__main__":
    main()
