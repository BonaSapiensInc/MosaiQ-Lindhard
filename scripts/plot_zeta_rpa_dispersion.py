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

"""Figure X candidate: Zeta-RPA plasmon rescue vs scalar RPA and Bohm–Gross.

Reads ``output/output_zeta_rpa_dispersion.dat`` from:

    ./simulator/build/mosaiq_simulator --pathway zeta-rpa --gamma GAMMA r_s T_kelvin

Columns: q  omega_p_RPA  omega_p_zeta  ImEps_RPA  ImEps_zeta  bohm_gross_estimate

Two-panel layout (mirrors ``plot_plasmon_dispersion.py``):
  (a) Dispersion: Bohm–Gross, scalar RPA, Zeta-RPA
  (b) Landau damping Im[ε] at each pathway's extracted pole
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any, cast

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import numpy as np

from plot_common import OUTPUT_DIR, apply_pdf_rcparams, format_temperature_k, save_figure

DEFAULT_DAT = OUTPUT_DIR / "output_zeta_rpa_dispersion.dat"
DEFAULT_FIG_STEM = "zeta_rpa_dispersion"

RS_HEADER = re.compile(r"r_s\s*=\s*([^\s]+)")
T_HEADER = re.compile(r"T_K\s*=\s*([^\s]+)")
GAMMA_HEADER = re.compile(r"Gamma_plasma\s*=\s*([^\s]+)")
PATHWAY_HEADER = re.compile(r"ResponsePathway\s*=\s*(\S+)")

COLOR_BOHM_GROSS = "#555555"
COLOR_RPA = "#8b1e1e"
COLOR_ZETA = "#0b5ea8"
COLOR_RPA_GAP = "#c47a7a"

PANEL_Q_SNAP = 0.5
PANEL_OMEGA_SNAP = 0.5
PANEL_Q_MAX_CAP = 4.0  # diagnostic focus; full .dat may extend farther


def configure_matplotlib() -> None:
    apply_pdf_rcparams(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 0.8,
            "grid.alpha": 0.30,
            "grid.linewidth": 0.5,
        }
    )


def finite_segments(q: np.ndarray, y: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split (q, y) into contiguous segments where y is finite."""
    mask = np.isfinite(y)
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    start: int | None = None

    for index, is_finite in enumerate(mask):
        if is_finite and start is None:
            start = index
        elif not is_finite and start is not None:
            segments.append((q[start:index], y[start:index]))
            start = None

    if start is not None:
        segments.append((q[start:], y[start:]))

    return segments


def snap_ceiling(value: float, step: float) -> float:
    return float(np.ceil(value / step) * step)


def parse_run_metadata(path: Path) -> tuple[float | None, float | None, float | None, str | None]:
    rs_value: float | None = None
    t_value: float | None = None
    gamma_value: float | None = None
    pathway: str | None = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            rs_match = RS_HEADER.search(line)
            if rs_match:
                rs_value = float(rs_match.group(1))
            t_match = T_HEADER.search(line)
            if t_match:
                t_value = float(t_match.group(1))
            gamma_match = GAMMA_HEADER.search(line)
            if gamma_match:
                gamma_value = float(gamma_match.group(1))
            pathway_match = PATHWAY_HEADER.search(line)
            if pathway_match:
                pathway = pathway_match.group(1)
    return rs_value, t_value, gamma_value, pathway


def load_dispersion(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if not path.is_file():
        raise SystemExit(
            f"Error: data file not found: {path}\n"
            "  Regenerate with:\n"
            "    ./simulator/build/mosaiq_simulator --pathway zeta-rpa --gamma 100 2.0 10000"
        )

    table = np.loadtxt(path, comments="#")
    if table.ndim != 2 or table.shape[1] < 6:
        raise SystemExit(f"Error: expected ≥6 columns in {path}, got shape {table.shape}")

    q = table[:, 0]
    omega_rpa = table[:, 1]
    omega_zeta = table[:, 2]
    im_rpa = table[:, 3]
    im_zeta = table[:, 4]
    bohm_gross = table[:, 5]
    return q, omega_rpa, omega_zeta, im_rpa, im_zeta, bohm_gross


def compute_display_bounds(
    q: np.ndarray,
    omega_rpa: np.ndarray,
    omega_zeta: np.ndarray,
    bohm_gross: np.ndarray,
    q_max_cap: float = PANEL_Q_MAX_CAP,
) -> tuple[float, float, float, float]:
    active = (
        (q > 0.0)
        & (q <= q_max_cap)
        & (np.isfinite(omega_rpa) | np.isfinite(omega_zeta) | np.isfinite(bohm_gross))
    )
    if not np.any(active):
        return 0.0, min(2.0, q_max_cap), 0.0, 4.0

    q_hi = min(q_max_cap, snap_ceiling(float(np.max(q[active])), PANEL_Q_SNAP))
    omega_stack = np.concatenate(
        [
            omega_rpa[active][np.isfinite(omega_rpa[active])],
            omega_zeta[active][np.isfinite(omega_zeta[active])],
            bohm_gross[active][np.isfinite(bohm_gross[active])],
        ]
    )
    if omega_stack.size == 0:
        return 0.0, q_hi, 0.0, 4.0

    w_hi = max(2.0, snap_ceiling(float(np.max(omega_stack)) * 1.08, PANEL_OMEGA_SNAP))
    return 0.0, q_hi, 0.0, w_hi


def mark_rpa_breakdown(
    ax: Axes,
    q: np.ndarray,
    omega_rpa: np.ndarray,
    omega_zeta: np.ndarray,
) -> None:
    """Highlight q where Zeta finds a pole but scalar RPA does not (rescue markers)."""
    rescued = np.isfinite(omega_zeta) & ~np.isfinite(omega_rpa)
    if not np.any(rescued):
        return
    ax.scatter(
        q[rescued],
        omega_zeta[rescued],
        s=28,
        marker="o",
        facecolors="none",
        edgecolors=COLOR_ZETA,
        linewidths=1.4,
        zorder=6,
        label="_nolegend_",
    )


def plot_rescue(
    q: np.ndarray,
    omega_rpa: np.ndarray,
    omega_zeta: np.ndarray,
    im_rpa: np.ndarray,
    im_zeta: np.ndarray,
    bohm_gross: np.ndarray,
    display_bounds: tuple[float, float, float, float],
    *,
    rs: float | None,
    t_kelvin: float | None,
    gamma: float | None,
    pathway: str | None,
) -> Figure:
    configure_matplotlib()

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(7.0, 5.8),
        gridspec_kw={"height_ratios": [1.45, 1.0], "hspace": 0.08},
    )

    q_lo, q_hi, w_lo, w_hi = display_bounds
    q_mask = (q >= q_lo) & (q <= q_hi)
    q_view = q[q_mask]

    for q_seg, y_seg in finite_segments(q_view, bohm_gross[q_mask]):
        ax_top.plot(
            q_seg,
            y_seg,
            color=COLOR_BOHM_GROSS,
            linewidth=1.8,
            linestyle="--",
            solid_capstyle="butt",
            dash_capstyle="butt",
            zorder=3,
        )

    for q_seg, y_seg in finite_segments(q_view, omega_rpa[q_mask]):
        ax_top.plot(
            q_seg,
            y_seg,
            color=COLOR_RPA,
            linewidth=2.0,
            linestyle="-",
            solid_capstyle="round",
            zorder=4,
        )

    for q_seg, y_seg in finite_segments(q_view, omega_zeta[q_mask]):
        ax_top.plot(
            q_seg,
            y_seg,
            color=COLOR_ZETA,
            linewidth=2.2,
            linestyle="-",
            solid_capstyle="round",
            zorder=5,
        )

    mark_rpa_breakdown(ax_top, q_view, omega_rpa[q_mask], omega_zeta[q_mask])

    # Soft band where RPA roots vanish but Zeta persists (visual “breakdown” cue).
    rpa_ok = np.isfinite(omega_rpa[q_mask])
    zeta_ok = np.isfinite(omega_zeta[q_mask])
    if np.any(rpa_ok) and np.any(zeta_ok & ~rpa_ok):
        q_last_rpa = float(np.max(q_view[rpa_ok]))
        ax_top.axvspan(
            q_last_rpa,
            q_hi,
            color=COLOR_RPA_GAP,
            alpha=0.12,
            zorder=1,
            linewidth=0,
        )
        ax_top.axvline(
            q_last_rpa,
            color=COLOR_RPA,
            linestyle=":",
            linewidth=1.2,
            alpha=0.85,
            zorder=2,
        )

    ax_top.set_ylabel(r"$\bar{\omega}_p = \hbar\omega_p/\epsilon_\mathrm{F}$")
    ax_top.set_xlim(q_lo, q_hi)
    ax_top.set_ylim(w_lo, w_hi)
    ax_top.margins(x=0, y=0)
    ax_top.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.35)

    legend_handles = [
        Line2D(
            [0], [0], color=COLOR_BOHM_GROSS, linewidth=1.8, linestyle="--"
        ),
        Line2D([0], [0], color=COLOR_RPA, linewidth=2.0),
        Line2D([0], [0], color=COLOR_ZETA, linewidth=2.2),
    ]
    legend_labels = [
        "Bohm--Gross",
        r"scalar RPA $\mathrm{Re}\,\varepsilon=0$",
        r"Zeta-RPA $\mathrm{Re}\,\varepsilon^\zeta=0$",
    ]
    if np.any(np.isfinite(omega_zeta[q_mask]) & ~np.isfinite(omega_rpa[q_mask])):
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markeredgecolor=COLOR_ZETA,
                markerfacecolor="none",
                markersize=7,
                linestyle="None",
                markeredgewidth=1.4,
            )
        )
        legend_labels.append("Zeta rescue (no RPA root)")

    legend = ax_top.legend(
        legend_handles,
        legend_labels,
        loc="upper left",
        frameon=True,
        framealpha=0.94,
        edgecolor="0.75",
        facecolor="white",
    )

    meta_bits: list[str] = []
    if rs is not None:
        meta_bits.append(rf"$r_s = {rs:g}$")
    if gamma is not None:
        meta_bits.append(rf"$\Gamma = {gamma:g}$")
    if t_kelvin is not None:
        meta_bits.append(format_temperature_k(t_kelvin))
    if pathway:
        meta_bits.append(pathway)

    fig.canvas.draw()
    canvas = cast(Any, fig.canvas)
    legend_bbox = legend.get_window_extent(canvas.get_renderer()).transformed(
        ax_top.transAxes.inverted()
    )
    if meta_bits:
        ax_top.text(
            legend_bbox.x0,
            legend_bbox.y0 - 0.014,
            r"$\;|\;$".join(meta_bits),
            transform=ax_top.transAxes,
            fontsize=9,
            va="top",
            ha="left",
            zorder=10,
            color="0.15",
        )

    # --- Panel (b): Landau damping at poles ---
    im_rpa_plot = np.where(np.isfinite(im_rpa) & (np.abs(im_rpa) > 0.0), np.abs(im_rpa), np.nan)
    im_zeta_plot = np.where(
        np.isfinite(im_zeta) & (np.abs(im_zeta) > 0.0), np.abs(im_zeta), np.nan
    )

    for q_seg, y_seg in finite_segments(q_view, im_rpa_plot[q_mask]):
        ax_bottom.plot(
            q_seg,
            y_seg,
            color=COLOR_RPA,
            linewidth=2.0,
            solid_capstyle="round",
            zorder=3,
        )

    for q_seg, y_seg in finite_segments(q_view, im_zeta_plot[q_mask]):
        ax_bottom.plot(
            q_seg,
            y_seg,
            color=COLOR_ZETA,
            linewidth=2.2,
            solid_capstyle="round",
            zorder=4,
        )

    ax_bottom.set_yscale("log")
    ax_bottom.set_ylabel(r"$|\Im\varepsilon(\bar{q}, \bar{\omega}_p)|$")
    ax_bottom.set_xlabel(r"Reduced wavevector $\bar{q} = q/k_\mathrm{F}$")
    ax_bottom.set_xlim(q_lo, q_hi)
    ax_bottom.margins(x=0)
    ax_bottom.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.35)

    ax_bottom.legend(
        [
            Line2D([0], [0], color=COLOR_RPA, linewidth=2.0),
            Line2D([0], [0], color=COLOR_ZETA, linewidth=2.2),
        ],
        [
            r"$|\Im\varepsilon^{\mathrm{RPA}}|$ at RPA pole",
            r"$|\Im\varepsilon^\zeta|$ at Zeta pole",
        ],
        loc="upper right",
        frameon=True,
        framealpha=0.94,
        edgecolor="0.75",
        facecolor="white",
        fontsize=9,
    )

    ax_top.text(
        0.98,
        0.98,
        "(a)",
        transform=ax_top.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="right",
        zorder=10,
    )
    ax_bottom.text(
        0.98,
        0.08,
        "(b)",
        transform=ax_bottom.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="right",
        zorder=10,
    )

    fig.align_ylabels([ax_top, ax_bottom])
    fig.subplots_adjust(hspace=0.08)
    return fig


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plot Zeta-RPA vs scalar RPA plasmon rescue (Figure 10)."
    )
    parser.add_argument(
        "--dat",
        type=Path,
        default=DEFAULT_DAT,
        help=f"dispersion comparison table (default: {DEFAULT_DAT})",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="optional output PDF path (default: output/zeta_rpa_dispersion.pdf + manuscript/figures/)",
    )
    parser.add_argument(
        "--q-max",
        type=float,
        default=PANEL_Q_MAX_CAP,
        help=f"display q ceiling (default: {PANEL_Q_MAX_CAP})",
    )
    args = parser.parse_args(argv)

    q, omega_rpa, omega_zeta, im_rpa, im_zeta, bohm_gross = load_dispersion(args.dat)
    rs, t_kelvin, gamma, pathway = parse_run_metadata(args.dat)
    bounds = compute_display_bounds(q, omega_rpa, omega_zeta, bohm_gross, args.q_max)

    fig = plot_rescue(
        q,
        omega_rpa,
        omega_zeta,
        im_rpa,
        im_zeta,
        bohm_gross,
        bounds,
        rs=rs,
        t_kelvin=t_kelvin,
        gamma=gamma,
        pathway=pathway,
    )

    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, format="pdf", dpi=300, bbox_inches="tight")
        saved_path = args.out
        plt.close(fig)
    else:
        saved_path = save_figure(fig, DEFAULT_FIG_STEM)
        plt.close(fig)

    n_rpa = int(np.sum(np.isfinite(omega_rpa)))
    n_zeta = int(np.sum(np.isfinite(omega_zeta)))
    n_rescue = int(np.sum(np.isfinite(omega_zeta) & ~np.isfinite(omega_rpa)))
    print(f"Loaded {q.size} q points from {args.dat}")
    print(f"  r_s={rs}  T={t_kelvin} K  Gamma={gamma}  pathway={pathway}")
    print(f"  roots: RPA={n_rpa}  Zeta={n_zeta}  Zeta-only rescue={n_rescue}")
    print(f"Saved {saved_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
