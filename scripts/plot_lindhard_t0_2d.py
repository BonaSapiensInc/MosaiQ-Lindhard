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
from matplotlib import cm

from plot_common import apply_pdf_rcparams, save_figure


def calc_t0_lindhard(q_grid: np.ndarray, w_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calculate analytic T=0 Lindhard function normalized by hbar * D(E_F)."""
    q = q_grid
    w = w_grid

    xi_minus = 0.5 * (w / q - q)
    xi_plus = 0.5 * (w / q + q)

    xi_minus_sq = xi_minus**2
    xi_plus_sq = xi_plus**2

    # --- Imaginary Part (Piecewise) ---
    im_chi = np.zeros_like(q)

    mask1 = (xi_minus_sq <= 1.0) & (xi_plus_sq >= 1.0)
    mask2 = (xi_minus_sq <= 1.0) & (xi_plus_sq <= 1.0)

    im_chi[mask1] = 1.0 - xi_minus_sq[mask1]
    im_chi[mask2] = w[mask2]
    im_chi *= np.pi / (4.0 * q)

    # --- Real Part (Logarithmic) ---
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
    cmap: cm.Colormap = cm.rainbow,
) -> None:
    apply_pdf_rcparams({"font.size": 12, "font.family": "serif"})
    fig, ax = plt.subplots(figsize=(8, 6))

    levels = np.linspace(np.min(z_grid), np.max(z_grid), 150)
    contour = ax.contourf(q_grid, w_grid, z_grid, levels=levels, cmap=cmap)

    line_levels = np.linspace(np.min(z_grid), np.max(z_grid), 25)
    ax.contour(
        q_grid,
        w_grid,
        z_grid,
        levels=line_levels,
        colors="white",
        linewidths=0.45,
        alpha=0.75,
    )

    cbar = fig.colorbar(contour, ax=ax, pad=0.02)
    cbar.set_label(label)

    ax.set_xlabel(r"Wave vector $\bar{q} = q/k_F$")
    ax.set_ylabel(r"Frequency $\bar{\omega} = \hbar\omega/\epsilon_F$")
    ax.set_title(title, fontweight="bold")

    ax.axvline(x=1.0, color="white", linestyle="--", linewidth=1.5, alpha=0.9)

    plt.tight_layout()
    saved_path = save_figure(fig, filename)
    plt.close(fig)
    print(f"Saved {saved_path}")


def main() -> None:
    q_1d = np.linspace(0.01, 4.0, 400)
    w_1d = np.linspace(0.01, 3.0, 400)
    q_grid, w_grid = np.meshgrid(q_1d, w_1d)

    re_chi, im_chi = calc_t0_lindhard(q_grid, w_grid)

    render_contour(
        q_grid,
        w_grid,
        re_chi,
        title=r"Analytic $T=0$ Lindhard: Real Part",
        label=r"$\Re\tilde{\chi}_{T=0}$ [$D(\epsilon_F)$ units]",
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
