# MosaiQ-Lindhard: Deterministic Linear Response Representation of Two-Fermion Plasmas

[![License: Custom Source-Available](https://img.shields.io/badge/License-Custom%20Source--Available-blue.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](simulator/)
[![CMake](https://img.shields.io/badge/Build-CMake-064F8C.svg)](simulator/CMakeLists.txt)
[![Python 3](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](scripts/)

**MosaiQ-Lindhard** is a C++20 engine and companion visualization toolchain for finite-temperature linear response in homogeneous two-component (electron–ion) plasmas. The engine evaluates Lindhard susceptibilities and two-component RPA dielectric functions **exactly on the real frequency axis** via Kramers–Kronig (Hilbert) relations and sinc-quadrature, bypassing divergent Matsubara pole summations entirely.

Companion manuscript: *Linear response representation of two-fermion plasmas* ([`manuscript/two-fermi.tex`](manuscript/two-fermi.tex)).

---

## Overview

Standard contour-integration and Matsubara formulations of finite-temperature Lindhard response accumulate systematic errors on the real axis—most visibly as spurious $\mathcal{O}(1)$ offsets in the real part of $\chi^0(q,\omega)$ at low temperature. MosaiQ-Lindhard instead constructs the causal response from its imaginary part through a deterministic Hilbert transform, with quadrature weights derived from sinc interpolation on a uniform grid. The result is a numerically stable, audit-ready linear-response engine suitable for structure-factor construction, sum-rule verification, and macroscopic ($q \lesssim 50$) plasmon spectroscopy.

---

## Key Features

- **Exact causality constraint** — Real and imaginary parts of $\chi^0_{ee}$, $\chi^0_{ii}$, and $\chi^0_{ei}$ are linked by the Hilbert transform; $\varepsilon(q,\omega)$ is assembled in the two-component RPA dielectric formalism without analytic continuation artifacts.
- **Elimination of $\mathcal{O}(1)$ real-part errors** — Removes the low-$T$ real-axis bias endemic to standard Matsubara/contour schemes; validated against sum rules and limiting cases (see `simulator/tests/`).
- **Dynamic structure factor contours** — $S_{ee}$, $S_{ii}$, $S_{ei}(\bar{q},\bar{\omega})$ logarithmic contour panels (Figure 3; $S_{ii}$ to $\bar{q}=50$).
- **Deterministic plasmon root-finding** — Brent-based pole extraction and $\gamma$-sweep static structure factors $S(q)$ up to $q = 50.0$.

---

## Repository Layout

| Path | Purpose |
|------|---------|
| [`simulator/`](simulator/) | C++20 engine (`mosaiq_simulator` CLI) and CTest suite |
| [`scripts/`](scripts/) | Python figure generators, figure sync, and manuscript packaging |
| [`manuscript/`](manuscript/) | LaTeX manuscript source (`two-fermi.tex`, `two-fermi.bib`; cover letter kept locally, not synced) |
| [`manuscript/figures/`](manuscript/figures/) | Manuscript figure PDFs (synced from `output/`) |
| [`manuscript/submission.zip`](manuscript/submission.zip) | Optional manuscript zip (`two-fermi.tex` + cited PDFs only; generated locally, not in git) |
| [`docs/`](docs/) | Architecture map ([`docs/README.md`](docs/README.md)); dual-pathway theory ([`theory.md`](docs/theory.md)), CLI ([`usage.md`](docs/usage.md)), Phase~2 Berezinian ([`BEREZINIAN_ARCHITECTURE.md`](docs/BEREZINIAN_ARCHITECTURE.md)) |
| [`output/`](output/) | Simulator grids (`.dat`, local only) and manuscript figures (`*.pdf`, tracked) |
| [`LICENSE`](LICENSE) | MosaiQ-Lindhard Source Code License Agreement |

### Scripts

| Script | Purpose |
|--------|---------|
| [`plot_lindhard_base.py`](scripts/plot_lindhard_base.py) | Figure 3 — bare Lindhard contours |
| [`plot_Sqw.py`](scripts/plot_Sqw.py) | Figure 4 — dynamic structure factor contours |
| [`plot_plasmon_dispersion.py`](scripts/plot_plasmon_dispersion.py) | Figure 5 — plasmon dispersion |
| [`plot_zeta_rpa_dispersion.py`](scripts/plot_zeta_rpa_dispersion.py) | Figure 10 — Zeta-RPA plasmon rescue |
| [`plot_lindhard_t0_2d.py`](scripts/plot_lindhard_t0_2d.py) | Figure 6 — $T=0$ contour panels |
| [`plot_t0_error.py`](scripts/plot_t0_error.py) | Figure 6 — $T=0$ validation slice |
| [`plot_thermal_anatomy.py`](scripts/plot_thermal_anatomy.py) | Figure 7 — thermal sharpening |
| [`plot_Sq_gamma.py`](scripts/plot_Sq_gamma.py) | Figure 8 — $\Gamma$ sweep static $S(q)$ |
| [`regenerate_figures_3_4.sh`](scripts/regenerate_figures_3_4.sh) | Batch regenerate Figs. 3–4 |
| [`regenerate_sq_gamma.sh`](scripts/regenerate_sq_gamma.sh) | Batch regenerate Fig. 8 data |
| [`sync_manuscript_figures.py`](scripts/sync_manuscript_figures.py) | Copy cited PDFs from `output/` to `manuscript/figures/` |
| [`package_aps_submission.sh`](scripts/package_aps_submission.sh) | Build `manuscript/submission.zip` (`two-fermi.tex` + cited PDFs) |

---

## Build & Run

### Prerequisites

- C++20 compiler (GCC 11+, Clang 14+, or Apple Clang 15+)
- CMake 3.20+
- Python 3.10+ with dependencies in [`requirements.txt`](requirements.txt) (for figures only)

### Build the simulator

```bash
cmake -S simulator -B simulator/build -DCMAKE_BUILD_TYPE=Release
cmake --build simulator/build
```

**Editor / clangd:** After the first CMake configure, optionally link the compilation database at the repo root so clangd picks up per-file build flags for `.cpp` sources:

```bash
ln -sf simulator/build/compile_commands.json compile_commands.json
```

Workspace settings already point clangd at `simulator/build` via `--compile-commands-dir`. Headers under `simulator/src/{core,physics,engine}/` and `simulator/tests/` each carry a committed `compile_flags.txt` (`-std=c++20`, `-I` to `src/`, `-Wno-pragma-once-outside-header`) so opening a header as a main file still resolves includes before CMake has been run. **When adding a new source subdirectory, copy an existing `compile_flags.txt` into it** — otherwise the Problems tab will report missing headers and false C++14/C++17 diagnostics.

**Do not enable `mitaki28.vscode-clang` in this workspace.** That extension runs bare `clang` without `-std=c++20` and floods Problems with false `std::numbers` / “C++17 extension” errors on headers such as `physics/Constants.hpp`. C++ diagnostics are owned exclusively by **clangd** (Microsoft/`anysphere` C/C++ IntelliSense is disabled in `.vscode/settings.json`).

Run the CLI from the repository root:

```bash
./simulator/build/mosaiq_simulator
```

### Run tests

```bash
ctest --test-dir simulator/build --output-on-failure
```

### Regenerate figures

**Version control:** Manuscript figure PDFs under `output/` and `manuscript/figures/` (for example `chi_*.pdf`, `S_*_contour.pdf`, `thermal_anatomy_*.pdf`) are tracked in git. Plot scripts write to `output/` and mirror each manuscript figure into `manuscript/figures/` automatically. To refresh the submission bundle without re-rendering, run `python3 scripts/sync_manuscript_figures.py`. **`output/*.dat` simulator grids are not** — they are listed in [`.gitignore`](.gitignore) and must be generated locally after clone. Plotting scripts call `mosaiq_simulator` automatically when a required `.dat` file is missing.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Build simulator (once)
cmake -S simulator -B simulator/build -DCMAKE_BUILD_TYPE=Release
cmake --build simulator/build --target mosaiq_simulator

# Figures 3 & 4 (bare Lindhard + RPA structure factors; writes output/*.dat then PDFs)
./scripts/regenerate_figures_3_4.sh 1 10000 2 10000

# Figure 7 — thermal anatomy at r_s = 1.0 (T = 70, 1000, 10000 K; panel (b) cold rung = 70 K)
python3 scripts/plot_thermal_anatomy.py

# Or plot only, if .dat files already exist locally:
python3 scripts/plot_lindhard_base.py
python3 scripts/plot_Sqw.py
python3 scripts/plot_plasmon_dispersion.py

# Static S(q) gamma sweep (Figure 8)
./scripts/regenerate_sq_gamma.sh 1 10,50,100,150
python3 scripts/plot_Sq_gamma.py

# T = 0 validation panels (Figure 6)
python3 scripts/plot_lindhard_t0_2d.py
python3 scripts/plot_t0_error.py   # requires output/output_t0_limit.dat from CTest

# Sync manuscript/figures/ (optional zip via package script)
python3 scripts/sync_manuscript_figures.py
./scripts/package_aps_submission.sh   # → manuscript/submission.zip
```

Typical single-run CLI (writes `output/output_lindhard_base.dat`, `output/output_structure_factor.dat`, and related files):

```bash
./simulator/build/mosaiq_simulator 1.0 10000
```

Or run the gamma sweep directly:

```bash
./simulator/build/mosaiq_simulator --gamma-sweep 1 10,50,100,150
python3 scripts/plot_Sq_gamma.py
```

### Compile the manuscript

From `manuscript/`, after all referenced `figures/*.pdf` exist (run `python3 scripts/sync_manuscript_figures.py` if needed):

```bash
cd manuscript && pdflatex two-fermi.tex && pdflatex two-fermi.tex
```

---

## License & Commercial Use

> **Important:** This repository is **not** released under MIT, Apache, or other permissive open-source licenses.

MosaiQ-Lindhard is distributed under the **[MosaiQ-Lindhard Source Code License Agreement](LICENSE)** (Copyright © 2026 Bona Sapiens, Inc.).

| Use case | Terms |
|----------|-------|
| **Personal & academic research** | **Free.** You may use, modify, and distribute the software for non-commercial personal or academic research. Publications must cite the accompanying manuscript (see [Citation](#citation)). |
| **Commercial entities** | **Paid license required.** Corporations, for-profit firms, and commercial R&D organizations must obtain a written commercial license before any use. |
| **Government & public institutions** | **Paid license required.** This includes national laboratories (e.g., LLNL, LANL, ORNL), military agencies, state/federal research offices, and publicly funded institutions deploying the engine operationally. |

For licensing inquiries, fee schedules, and contract negotiation:

**Bona Sapiens, Inc.** — [kim.ingee@bonasapiens.com](mailto:kim.ingee@bonasapiens.com)

---

## Citation

If you use MosaiQ-Lindhard in academic work—including the causality-constrained Lindhard pipeline, Zeta-RPA, or the $DX=C$ Berezinian modules—please cite the permanent Zenodo archive:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21465936.svg)](https://doi.org/10.5281/zenodo.21465936)

> In-Gee Kim (2026). *Linear response representation of two-fermion plasmas & MosaiQ-Lindhard Engine* (v3.0.0-berezinian-freeze). Zenodo. https://doi.org/10.5281/zenodo.21465936

Manuscript source: [`manuscript/two-fermi.tex`](manuscript/two-fermi.tex).

---

## Data Availability

Source code, test suites, and build instructions are published under the **MosaiQ-Lindhard** name to avoid confusion with unrelated third-party packages.

**Repository:** [https://github.com/BonaSapiensInc/MosaiQ-Lindhard](https://github.com/BonaSapiensInc/MosaiQ-Lindhard)

**Zenodo archive:** [https://doi.org/10.5281/zenodo.21465936](https://doi.org/10.5281/zenodo.21465936) (release `v3.0.0-berezinian-freeze`)

**What is in git:** C++ engine, Python plotting scripts, LaTeX manuscript, and regenerated **figure PDFs** under `output/`.

**What is not in git:** **`output/*.dat`** numerical grids from the CLI (ignored by design) and **`manuscript/submission.zip`** (regenerate with [`package_aps_submission.sh`](scripts/package_aps_submission.sh)). After cloning, build `mosaiq_simulator` and run the [Regenerate figures](#regenerate-figures) commands above; each workflow documents the `(r_s, T)` or `--gamma-sweep` parameters used in the manuscript.
