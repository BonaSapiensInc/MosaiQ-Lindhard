#!/usr/bin/env python3
"""Compare scalar RPA vs Zeta-RPA (production + experimental) along q at fixed omega.

Reads output/zeta_rpa_scalar_compare.dat from verify_zeta_rpa_scalar.
Writes output/zeta_rpa_scalar_compare.pdf (and mirrors to manuscript/figures/ if present).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DAT = ROOT / "output" / "zeta_rpa_scalar_compare.dat"
OUT = ROOT / "output" / "zeta_rpa_scalar_compare.pdf"


def main() -> int:
    if not DAT.is_file():
        print(f"Missing {DAT}; run verify_zeta_rpa_scalar first.", file=sys.stderr)
        return 1

    data = np.loadtxt(DAT, comments="#")
    # columns: q omega ReChiL ImChiL ReChiRPA ImChiRPA ReChiZeta ImChiZeta ReChiExp ImChiExp Wprod Wexp
    omega0 = 0.5
    sl = np.isclose(data[:, 1], omega0)
    if not np.any(sl):
        omega0 = float(data[0, 1])
        sl = np.isclose(data[:, 1], omega0)

    q = data[sl, 0]
    rpa = data[sl, 4] + 1j * data[sl, 5]
    zeta = data[sl, 6] + 1j * data[sl, 7]
    exp = data[sl, 8] + 1j * data[sl, 9]

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), constrained_layout=True)
    axes[0].plot(q, np.real(rpa), "k-", lw=1.8, label=r"scalar RPA")
    axes[0].plot(q, np.real(zeta), "C0--", lw=1.6, label=r"ZetaRPA ($W=1$)")
    axes[0].plot(q, np.real(exp), "C3-", lw=1.4, label=r"ZetaRPA\_Experimental")
    axes[0].set_xlabel(r"$\bar{q}$")
    axes[0].set_ylabel(r"$\mathrm{Re}\,\chi$")
    axes[0].set_title(rf"$\bar{{\omega}}={omega0:g}$")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(q, np.imag(rpa), "k-", lw=1.8, label=r"scalar RPA")
    axes[1].plot(q, np.imag(zeta), "C0--", lw=1.6, label=r"ZetaRPA ($W=1$)")
    axes[1].plot(q, np.imag(exp), "C3-", lw=1.4, label=r"ZetaRPA\_Experimental")
    axes[1].set_xlabel(r"$\bar{q}$")
    axes[1].set_ylabel(r"$\mathrm{Im}\,\chi$")
    axes[1].set_title(rf"$\bar{{\omega}}={omega0:g}$")
    axes[1].legend(frameon=False, fontsize=8)

    fig.suptitle(r"Phase Z1 scalar: RPA vs Zeta-RPA ($\Gamma=50$, $r_s=1$)", fontsize=11)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
