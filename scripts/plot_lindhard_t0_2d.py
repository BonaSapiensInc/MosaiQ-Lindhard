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

"""Render T=0 analytic Lindhard function (Real and Imaginary) 2D contour maps."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize

from plot_common import save_figure
from plot_Sqw import (
    CONTOUR_ISOLINE_ALPHA,
    CONTOUR_ISOLINE_COLOR,
    CONTOUR_ISOLINE_LINEWIDTH,
    ELECTRON_DISPLAY_BOUNDS,
    configure_seaborn,
    filter_contour_line_levels,
    format_contour_display_axes,
    frame_contour_axis,
    log_contour_levels,
    prepare_plot_grid,
    rainbow_log_cmap,
)


def calc_t0_lindhard(q_grid: np.ndarray, w_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calculate analytic T=0 Lindhard function normalized by hbar * D(E_F)."""
    q = q_grid
    w = w_grid

    xi_minus = 0.5 * (w / q - q)
    xi_plus = 0.5 * (w / q + q)

    xi_minus_sq = xi_minus**2
    xi_plus_sq = xi_plus**2

    im_chi = np.zeros_like(q)

    mask1 = (xi_minus_sq <= 1.0) & (xi_plus_sq >= 1.0)
    mask2 = (xi_minus_sq <= 1.0) & (xi_plus_sq <= 1.0)

    im_chi[mask1] = 1.0 - xi_minus_sq[mask1]
    im_chi[mask2] = w[mask2]
    im_chi *= np.pi / (4.0 * q)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        term1 = q
        term2 = 0.5 * (1.0 - xi_minus_sq) * np.log(np.abs((xi_minus - 1.0) / (xi_minus + 1.0)))
        term3 = 0.5 * (1.0 - xi_plus_sq) * np.log(np.abs((xi_plus - 1.0) / (xi_plus + 1.0)))

        re_chi = (1.0 / (2.0 * q)) * (term1 + term2 + term3)
        re_chi = np.nan_to_num(re_chi, nan=0.0, posinf=0.0, neginf=0.0)

    return re_chi, im_chi


def render_contour(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    z_grid: np.ndarray,
    title: str,
    label: str,
    filename: str,
) -> None:
    """Match Figure 3 structure-factor contour styling (log rainbow + white isolines)."""
    configure_seaborn()

    plot_grid, norm = prepare_plot_grid(z_grid)
    color_norm: LogNorm | Normalize = norm if norm is not None else Normalize()
    cmap = rainbow_log_cmap()
    scalar_map = ScalarMappable(cmap=cmap, norm=color_norm)
    scalar_map.set_array([])

    fig, ax = plt.subplots(figsize=(10, 7), layout="constrained")
    if isinstance(color_norm, LogNorm):
        fill_levels, line_levels = log_contour_levels(color_norm, n_filled=200, n_lines=72)
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
        ax.contour(
            q_grid,
            w_grid,
            plot_grid,
            levels=line_levels,
            norm=color_norm,
            colors=CONTOUR_ISOLINE_COLOR,
            linewidths=CONTOUR_ISOLINE_LINEWIDTH,
            alpha=CONTOUR_ISOLINE_ALPHA,
        )
    else:
        ax.contourf(q_grid, w_grid, plot_grid, levels=150, cmap=cmap, norm=color_norm)
        ax.contour(
            q_grid,
            w_grid,
            plot_grid,
            levels=60,
            norm=color_norm,
            colors=CONTOUR_ISOLINE_COLOR,
            linewidths=CONTOUR_ISOLINE_LINEWIDTH,
            alpha=CONTOUR_ISOLINE_ALPHA,
        )

    frame_contour_axis(ax)
    ax.set_xlabel(r"Wave vector $q \ [k_F]$")
    ax.set_ylabel(r"Frequency $\omega \ [E_F/\hbar]$")
    ax.set_title(f"{title} — Contour Map", fontweight="bold")
    ax.axvline(x=1.0, color="white", linestyle="--", linewidth=1.5, alpha=0.9)
    format_contour_display_axes(
        ax,
        ELECTRON_DISPLAY_BOUNDS[0],
        ELECTRON_DISPLAY_BOUNDS[1],
        ELECTRON_DISPLAY_BOUNDS[2],
        ELECTRON_DISPLAY_BOUNDS[3],
    )
    fig.colorbar(scalar_map, ax=ax, fraction=0.046, pad=0.04, label=label)

    saved_path = save_figure(fig, filename)
    plt.close(fig)
    print(f"Saved {saved_path}")


def main() -> None:
    _, q_hi, _, w_hi = ELECTRON_DISPLAY_BOUNDS
    q_1d = np.linspace(0.01, q_hi, 400)
    w_1d = np.linspace(0.01, w_hi, 400)
    q_grid, w_grid = np.meshgrid(q_1d, w_1d)

    re_chi, im_chi = calc_t0_lindhard(q_grid, w_grid)

    render_contour(
        q_grid,
        w_grid,
        np.abs(re_chi),
        title=r"Analytic $T=0$ Lindhard: Real Part",
        label=r"$|\Re\tilde{\chi}_{T=0}|$ [$D(\epsilon_F)$ units]",
        filename="t0_analytic_re_contour",
    )
    render_contour(
        q_grid,
        w_grid,
        -im_chi,
        title=r"Analytic $T=0$ Lindhard: Imaginary Part (Dissipation)",
        label=r"$|\Im\tilde{\chi}_{T=0}|$ [$D(\epsilon_F)$ units]",
        filename="t0_analytic_im_contour",
    )


if __name__ == "__main__":
    main()
