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

"""Publication-quality T=0 Lindhard limit validation figure from MosaiQ-Lindhard output (PRE style)."""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np

from plot_common import OUTPUT_DIR, apply_pdf_rcparams, output_path, save_figure

DATA_FILE = OUTPUT_DIR / "output_t0_limit.dat"
OUTPUT_FIG = output_path("t0_limit_validation")

COLOR_RE = "#1f3b73"
COLOR_IM = "#8b1e1e"
COLOR_ERR_RE = "#1f3b73"
COLOR_ERR_IM = "#8b1e1e"

RS_HEADER = re.compile(r"r_s\s*=\s*([^\s]+)")
T_HEADER = re.compile(r"equivalent\s+T_K\s*~\s*([^\s]+)")


def load_t0_limit(path: Path) -> dict[str, np.ndarray]:
    if not path.is_file():
        raise SystemExit(f"Error: data file not found: {path}")

    table = np.loadtxt(path, comments="#")
    if table.ndim != 2 or table.shape[1] < 8:
        raise SystemExit(f"Error: expected 8 columns in {path}, got shape {table.shape}")

    return {
        "q": table[:, 0],
        "omega": table[:, 1],
        "re_exact": table[:, 2],
        "im_exact": table[:, 3],
        "re_engine": table[:, 4],
        "im_engine": table[:, 5],
        "error_re": table[:, 6],
        "error_im": table[:, 7],
    }


def parse_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if "q_bar =" in line:
                for token in line.replace("#", "").split():
                    if "=" in token:
                        key, value = token.split("=", 1)
                        metadata[key.strip()] = value.strip()
            if match := RS_HEADER.search(line):
                metadata["r_s"] = match.group(1)
            if match := T_HEADER.search(line):
                metadata["T_K"] = match.group(1)
    return metadata


def format_rs_label(r_s_raw: str | None) -> str:
    """Format r_s for figure titles (Fig. 6 uses r_s = 2.0 at tau = 10^{-4})."""
    if r_s_raw is None:
        return "2.0"
    return f"{float(r_s_raw):.1f}"


def format_temperature_label(t_k_raw: str | None) -> str:
    if t_k_raw is None:
        return r"11\,\mathrm{K}"
    return rf"{round(float(t_k_raw)):g}\,\mathrm{{K}}"


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
            "grid.alpha": 0.25,
            "grid.linewidth": 0.5,
        }
    )


def plot_t0_validation(
    data: dict[str, np.ndarray], rs_label: str, q_label: str, t_label: str
) -> Figure:
    configure_matplotlib()

    omega = data["omega"]

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(6.5, 5.5),
        gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.08},
    )

    ax_top.plot(
        omega,
        data["re_exact"],
        color=COLOR_RE,
        linewidth=1.8,
        label=r"$\Re\tilde{\chi}_{T=0}$ (exact)",
    )
    ax_top.plot(
        omega,
        data["im_exact"],
        color=COLOR_IM,
        linewidth=1.8,
        linestyle="--",
        dashes=(5, 3),
        label=r"$\Im\tilde{\chi}_{T=0}$ (exact)",
    )
    ax_top.scatter(
        omega,
        data["re_engine"],
        s=14,
        facecolors="none",
        edgecolors=COLOR_RE,
        linewidths=0.9,
        marker="o",
        label=r"Engine $\Re\tilde{\chi}$",
        zorder=3,
    )
    ax_top.scatter(
        omega,
        data["im_engine"],
        s=14,
        facecolors="none",
        edgecolors=COLOR_IM,
        linewidths=0.9,
        marker="s",
        label=r"Engine $\Im\tilde{\chi}$",
        zorder=3,
    )

    ax_top.set_ylabel(r"$\tilde{\chi}(\bar{q}, \bar{\omega})$ [$D(\epsilon_\mathrm{F})$ units]")
    ax_top.set_title(
        rf"Highly degenerate limit ($r_s = {rs_label}$) at $\bar{{q}} = {q_label}$ at $T = {t_label}$",
        fontsize=11,
    )
    ax_top.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.35)
    ax_top.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="0.75")

    positive_re = np.maximum(data["error_re"], np.finfo(float).tiny)
    positive_im = np.maximum(data["error_im"], np.finfo(float).tiny)

    ax_bottom.plot(
        omega,
        positive_re,
        color=COLOR_ERR_RE,
        linewidth=1.6,
        label=r"$|\Delta\Re|$",
    )
    ax_bottom.plot(
        omega,
        positive_im,
        color=COLOR_ERR_IM,
        linewidth=1.6,
        linestyle="--",
        dashes=(5, 3),
        label=r"$|\Delta\Im|$",
    )

    ax_bottom.set_yscale("log")
    ax_bottom.set_ylabel(r"Absolute error")
    ax_bottom.set_xlabel(r"Reduced frequency $\bar{\omega} = \hbar\omega / \epsilon_\mathrm{F}$")
    ax_bottom.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.35)
    ax_bottom.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="0.75")

    fig.align_ylabels([ax_top, ax_bottom])
    fig.subplots_adjust(hspace=0.08)
    return fig


def main() -> None:
    metadata = parse_metadata(DATA_FILE)
    data = load_t0_limit(DATA_FILE)

    q_values = np.unique(data["q"])
    if q_values.size != 1:
        raise SystemExit(
            f"Error: expected a single q_bar in {DATA_FILE}, found {q_values.size} values."
        )

    rs_label = format_rs_label(metadata.get("r_s"))
    q_label = metadata.get("q_bar", f"{q_values[0]:g}")
    t_label = format_temperature_label(metadata.get("T_K"))
    fig = plot_t0_validation(data, rs_label, q_label, t_label)

    finite_re_errors = data["error_re"][np.isfinite(data["error_re"])]
    finite_im_errors = data["error_im"][np.isfinite(data["error_im"])]

    saved_path = save_figure(fig, OUTPUT_FIG.stem)
    plt.close(fig)

    print(f"Loaded {data['omega'].size} omega points at q_bar = {q_values[0]:g}.")
    if finite_re_errors.size:
        print(f"max |Delta Re| = {np.max(finite_re_errors):.3e}")
    else:
        print("max |Delta Re| = n/a (no finite reference points)")
    print(f"max |Delta Im| = {np.max(finite_im_errors):.3e}")
    print(f"Saved {saved_path}")


if __name__ == "__main__":
    main()
