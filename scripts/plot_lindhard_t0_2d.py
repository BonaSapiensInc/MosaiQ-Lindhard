#!/usr/bin/env python3
"""Render T=0 analytic Lindhard function (Real and Imaginary) 2D contour maps."""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"


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
    plt.rcParams.update({"font.size": 12, "font.family": "serif"})
    fig, ax = plt.subplots(figsize=(8, 6))

    levels = np.linspace(np.min(z_grid), np.max(z_grid), 150)
    contour = ax.contourf(q_grid, w_grid, z_grid, levels=levels, cmap=cmap)

    line_levels = np.linspace(np.min(z_grid), np.max(z_grid), 25)
    ax.contour(
        q_grid,
        w_grid,
        z_grid,
        levels=line_levels,
        colors="black",
        linewidths=0.8,
        alpha=0.8,
    )

    cbar = fig.colorbar(contour, ax=ax, pad=0.02)
    cbar.set_label(label)

    ax.set_xlabel(r"Wave vector $\bar{q} = q/k_F$")
    ax.set_ylabel(r"Frequency $\bar{\omega} = \hbar\omega/\epsilon_F$")
    ax.set_title(title, fontweight="bold")

    ax.axvline(x=1.0, color="white", linestyle="--", linewidth=1.5, alpha=0.9)
    ax.text(
        1.05,
        np.max(w_grid) * 0.9,
        r"$\bar{q}=1$ slice",
        color="white",
        fontsize=11,
        fontweight="bold",
    )

    plt.tight_layout()
    output_path = OUTPUT_DIR / filename
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}")


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
        filename="t0_analytic_re_contour.pdf",
    )
    render_contour(
        q_grid,
        w_grid,
        -im_chi,
        title=r"Analytic $T=0$ Lindhard: Imaginary Part (Dissipation)",
        label=r"$|\Im\tilde{\chi}_{T=0}|$ [$D(\epsilon_F)$ units]",
        filename="t0_analytic_im_contour.pdf",
    )


if __name__ == "__main__":
    main()
